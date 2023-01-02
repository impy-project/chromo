import multiprocessing as mp
import traceback
import sys

__all__ = ["make_remote_controlled_model"]


def run(output, Model, seed, init_kwargs, stables, random_state, method, args):
    model = Model(seed, **init_kwargs)
    if random_state is not None:
        model.random_state = random_state
    try:
        if method == "call":
            for k in stables - model.get_stable():
                model.set_stable(k, True)
            for k in model.get_stable() - stables:
                model.set_stable(k, False)
            for x in model(*args):
                output.put(x)
        else:
            x = getattr(model, method)(*args)
            output.put(x)
    except Exception as exc:
        (msg,) = exc.args
        tb = sys.exc_info()[2]
        s = "".join(traceback.format_tb(tb))
        exc.args = (f"{msg}\n\nBacktrace from child process:\n{s}",)
        output.put(exc)
        return
    output.put(model.random_state)


def get(timeout, process, output):
    from queue import Empty

    for _ in range(timeout):
        if not process.is_alive():
            raise ValueError("process died")
        try:
            x = output.get(timeout=1)
            break
        except Empty:
            pass
    else:
        raise TimeoutError("process send no data")

    if isinstance(x, Exception):
        raise x
    return x


def make_remote_controlled_model(name, Model):
    """
    Dynamically create a Model class which executes calls to Model.cross_section and
    Model.__call__ to the wrapped model in a separate process to work around
    initialization issues.

    Parameters
    ----------
    name : str
        Name of the generated class, must be unique.
    Model : class
        The wrapped class.
    """
    # This is a big hack, but required until a proper solution is available. A new process
    # is created whenever the cross_section or __call__ methods are used.
    #
    # Generating events is optimized, the process is kept alive during generation and
    # sends its events via a queue to the main process. In an exception occurs in the
    # child process, it is passed to the main process and raised there.
    #
    # Sending of data currently uses pickle, which has some overhead. This penalty only
    # matters for event generation, where the performance drop could be noticable. We
    # could avoid this by allocating EventData in shared memory, but putting more time and
    # effort into this for a model that is barely used (Phojet) does not seem effective.
    #
    # To maintain the illusion that we are pulling events from a model instance instead of
    # restarting it repeatedly, we roundtrip the state of the random number generators
    # everytime we create a new process.

    def __init__(self, seed=None, timeout=100, **kwargs):
        Model.__init__(self, seed, **kwargs)
        self._timeout = timeout
        self._init_kwargs = kwargs

    def __call__(self, kin, nevents):
        process, output = self._remote_init("call", (kin, nevents))
        for _ in range(nevents):
            yield get(self._timeout, process, output)
        self.random_state = get(self._timeout, process, output)
        process.join()

    def _cross_section(self, kin):
        return self._remote_method("_cross_section", kin)

    def _remote_init(self, method, args):
        ctx = mp.get_context("spawn")
        output = ctx.Queue()
        process = ctx.Process(
            target=run,
            args=(
                output,
                Model,
                self.seed,
                self._init_kwargs,
                self.get_stable(),
                self.random_state,
                method,
                args,
            ),
            daemon=True,
        )
        process.start()
        return process, output

    def _remote_method(self, method, *args):
        process, output = self._remote_init(method, args)
        x = get(self._timeout, process, output)
        self.random_state = get(self._timeout, process, output)
        process.join()
        return x

    # dynamically create a new class which inherits from Model and
    # overrides the methods which need to do remote calls
    return type(
        name,
        (Model,),
        {
            "__init__": __init__,
            "__call__": __call__,
            "_cross_section": _cross_section,
            "_remote_init": _remote_init,
            "_remote_method": _remote_method,
        },
    )
