import multiprocessing as mp
import traceback
import sys
from impy.common import MCRun, RMMARDState
from queue import Empty


__all__ = ["MCRunRemote"]


class MCRunRemote(MCRun):
    """
    Base class for models which cannot switch beam particles arbitrarily.

    Calls to Model.cross_section and Model.__call__ are executed
    in a separate process with a fresh instance of the underlying generator
    to work around initialization issues.

    If the parameter timeout is set to zero or negative values, the remote
    control feature is disabled. This may be required when you want to run
    the models in parallel in a Notebook with joblib, and possibly elsewhere.
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

    def __init__(self, seed=None, *, timeout=100, **kwargs):
        super().__init__(seed, timeout=timeout, **kwargs)

    def __call__(self, kin, nevents):
        if self._timeout > 0:
            with _RemoteCall(self, "call", kin, nevents) as rc:
                for _ in range(nevents):
                    yield rc.get()
        else:
            return super().__call__(kin, nevents)

    def cross_section(self, kin, **kwargs):
        if self._timeout > 0:
            with _RemoteCall(self, "cross_section", kin, **kwargs) as rc:
                return rc.get()
        else:
            return super().cross_section(kin, **kwargs)


class _RemoteCall:
    def __init__(self, parent, method, *args, **kwargs):
        self.parent = parent
        self.timeout = parent._timeout
        ctx = mp.get_context("spawn")
        # Queue should only hold one item at a time
        self.output = ctx.Queue(maxsize=1)
        self.stop = ctx.Event()
        state = (
            parent.seed,
            parent._init_kwargs,
            parent.get_stable(),
            parent.random_state,
        )
        self.process = ctx.Process(
            target=_run,
            args=(
                self.output,
                self.stop,
                self.parent.__class__,
                state,
                method,
                args,
                kwargs,
            ),
            daemon=True,
        )
        self.process.start()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, tb):
        # maybe stop iteration in remote process
        self.stop.set()
        # empty queue until we get RMMARDState
        while self.output.full():
            x = self.output.get()
            if isinstance(x, RMMARDState):
                self.parent.random_state = x
                break
        assert self.output.empty()
        self.output.close()
        self.process.join(timeout=1)
        if self.process.is_alive():
            self.process.kill()
        self.process.join()
        self.process.close()

    def get(self):
        for _ in range(self.timeout):
            if not self.process.is_alive():
                raise RuntimeError("process died")
            try:
                x = self.output.get(timeout=1)
                break
            except Empty:
                pass
        else:
            raise TimeoutError("process send no data")

        if isinstance(x, RemoteException):
            exc, stb = x.args
            (msg,) = exc.args
            exc.args = (f"{msg}\n\nBacktrace from child process:\n{stb}",)
            raise exc

        return x


class RemoteException(Exception):
    pass


def _run(output, stop, Model, state, method, args, kwargs):
    seed, init_kwargs, stable_state, random_state = state
    # timeout = 0 disables remote call, we don't want this to be recursive
    model = Model(seed, timeout=0, **init_kwargs)
    model.random_state = random_state
    try:
        if method == "call":
            for k in stable_state - model.get_stable():
                model.set_stable(k, True)
            for k in model.get_stable() - stable_state:
                model.set_stable(k, False)
            for x in model(*args, **kwargs):
                # If iteration is stopped in the main thread,
                # we must stop it also here, to get our
                # random_state below
                if stop.is_set():
                    break
                # The copy is necessary, but shouldn't be. This
                # smells like a bug in the pickling of MCEvent, or
                # maybe the process.Queue does not use pickling.
                output.put(x.copy())
        else:
            x = getattr(model, method)(*args, **kwargs)
            output.put(x)

    except Exception as exc:
        # To aid debugging, we catch the traceback in the
        # Remote process and pass it along with the exception
        # to the main process. Unfortunately, one cannot pickle
        # traceback objects (without external library support),
        # so we convert the traceback into a string here, and
        # append that string to the exception message later.
        tb = sys.exc_info()[2]
        s = "".join(traceback.format_tb(tb))
        exc2 = RemoteException(exc, s)
        output.put(exc2)

    # Finally also return random state, whether there
    # was an exception or not. If we switch to the Numpy
    # generator, we can probably allocate the random state
    # in shared memory, so it is shared at no additional cost,
    # and this manual synchronisation is no longer necessary.
    output.put(model.random_state)
