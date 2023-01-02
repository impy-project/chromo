import multiprocessing as mp
import traceback
import sys
from impy.common import MCRun

__all__ = ["MCRunRemote"]


def run(output, Model, seed, init_kwargs, stables, random_state, method, args):
    model = Model(seed, **init_kwargs)
    model.random_state = random_state
    try:
        if method == "call":
            for k in stables - model.get_stable():
                model.set_stable(k, True)
            for k in model.get_stable() - stables:
                model.set_stable(k, False)
            # we must call MCRun.__call__ here,
            # not MCRunRemote.__call__
            for x in MCRun.__call__(model, *args):
                output.put(x.copy())
        elif method == "cross_section":
            # same as above
            x = MCRun.cross_section(model, *args)
            output.put(x)
        else:
            assert False  # never arrive here
    except Exception as exc:
        (msg,) = exc.args
        tb = sys.exc_info()[2]
        s = "".join(traceback.format_tb(tb))
        exc.args = (f"{msg}\n\nBacktrace from child process:\n{s}",)
        output.put(exc)
        return
    # also return random state
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


class MCRunRemote(MCRun):
    """
    Base class for models which cannot switch beam particles arbitrarily.

    Calls to Model.cross_section and Model.__call__ are executed
    in a separate process with a fresh instance of the underlying generator
    to work around initialization issues.
    """

    # This is a big hack, but required until a proper solution is available.
    #
    # Generating events is optimized, the process is kept alive during generation and
    # sends its events via a queue to the main process. In an exception occurs in the
    # child process, it is passed to the main process and raised there.
    #
    # Sending of data may use pickle (seems to be platform dependent) and should have some
    # overhead. This penalty only matters for event generation, where the performance drop
    # could be noticable. We could avoid this by allocating EventData in shared memory.
    #
    # To maintain the illusion that we are pulling events from a model instance instead of
    # restarting it repeatedly, we must roundtrip the state of the random number generator
    # everytime we create a new process. So this functionality relies strongly on correct
    # RPNG state persistence.

    def __init__(self, seed=None, timeout=100, **kwargs):
        super().__init__(seed, **kwargs)
        self._timeout = timeout
        self._init_kwargs = kwargs

    def __call__(self, kin, nevents):
        process, output = self._remote_init("call", (kin, nevents))
        for _ in range(nevents):
            yield get(self._timeout, process, output)
        self.random_state = get(self._timeout, process, output)
        process.join()

    def cross_section(self, kin, **kwargs):
        return self._remote_method("cross_section", kin)

    def _remote_init(self, method, args):
        ctx = mp.get_context("spawn")
        output = ctx.Queue()
        process = ctx.Process(
            target=run,
            args=(
                output,
                self.__class__,
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
