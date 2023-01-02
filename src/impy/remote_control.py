import multiprocessing as mp
import traceback
import sys
from impy.common import MCRun
from queue import Empty


__all__ = ["MCRunRemote"]


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
        self._remote = _RemoteCall(self, timeout)

    def __call__(self, kin, nevents):
        with self._remote("call", kin, nevents) as rc:
            for _ in range(nevents):
                yield rc.get()

    def cross_section(self, kin, **kwargs):
        with self._remote("cross_section", kin, **kwargs) as rc:
            return rc.get()


class _RemoteCall:
    def __init__(self, parent, timeout):
        self.parent = parent
        self.timeout = timeout

    def __enter__(self):
        return self

    def __call__(self, method, *args, **kwargs):
        ctx = mp.get_context("spawn")
        self.output = ctx.Queue()
        p = self.parent
        state = (p.seed, p._init_kwargs, p.get_stable(), p.random_state)
        self.process = ctx.Process(
            target=_run,
            args=(
                self.output,
                self.parent.__class__,
                state,
                method,
                args,
                kwargs,
            ),
            daemon=True,
        )
        self.process.start()
        return self

    def __exit__(self, exc_type, exc_val, tb):
        if exc_type is None:
            self.parent.random_state = self.get()
        self.process.join()

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
            if self.process.is_alive():
                self.process.terminate()
            raise TimeoutError("process send no data")

        if isinstance(x, RemoteException):
            exc, stb = x.args
            (msg,) = exc.args
            exc.args = (f"{msg}\n\nBacktrace from child process:\n{stb}",)
            raise exc
        return x


class RemoteException(Exception):
    pass


def _run(output, Model, state, method, args, kwargs):
    seed, init_kwargs, stable_state, random_state = state
    model = Model(seed, **init_kwargs)
    model.random_state = random_state
    try:
        if method == "call":
            for k in stable_state - model.get_stable():
                model.set_stable(k, True)
            for k in model.get_stable() - stable_state:
                model.set_stable(k, False)
            # we must explicitly call MCRun.__call__ here,
            # to get MCRunRemote.__call__
            for x in MCRun.__call__(model, *args, **kwargs):
                output.put(x.copy())
        else:
            # same as above
            meth = getattr(MCRun, method)
            x = meth(model, *args, **kwargs)
            output.put(x)
    except Exception as exc:
        tb = sys.exc_info()[2]
        s = "".join(traceback.format_tb(tb))
        exc2 = RemoteException(exc, s)
        output.put(exc2)
    # also return random state
    output.put(model.random_state)
