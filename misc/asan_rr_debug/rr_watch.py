import gdb
import sys

def p(msg):
    print(msg)
    sys.stdout.flush()

class CrashBreak(gdb.Breakpoint):
    def stop(self):
        frame = gdb.selected_frame()
        try:
            thisval = frame.read_var("this")
        except Exception as e:
            p("could not read this: %s" % e)
            return True
        thisaddr = int(thisval)
        p("\n### CRASH SITE HIT: this=0x%x ###" % thisaddr)
        dataptr = None
        for expr in ["(long)this->d_X.data()", "(long)d_X.data()",
                     "(long)this->d_X.m_data", "(long)d_X.m_data"]:
            try:
                dataptr = gdb.parse_and_eval(expr)
                p("read via expr: %s" % expr)
                break
            except Exception as e:
                p("expr failed: %s %s" % (expr, e))
        if dataptr is None:
            p("could not read d_X data pointer via any expression")
            try:
                gdb.execute("print *this")
            except Exception:
                pass
            return True
        addr = int(dataptr)
        p("d_X data pointer: 0x%x" % addr)

        self.enabled = False  # avoid self-deletion issues; just disable
        gdb.post_event(lambda: walker.arm_and_go(addr))
        return True

class ReverseWalker(object):
    def __init__(self, n):
        self.remaining = n
        self.wp = None
        self.step = 0

    def arm_and_go(self, addr):
        self.addr = addr
        self.next_watch()

    def next_watch(self):
        if self.remaining <= 0:
            p("=== REVERSE WALK DONE (exhausted budget) ===")
            self.finish()
            return
        self.remaining -= 1
        self.step += 1
        wp_expr = "*(long*)0x%x" % self.addr
        try:
            gdb.execute("watch %s" % wp_expr)
            self.armed = True
        except Exception as e:
            p("watch failed: %s" % e)
            self.finish()
            return
        try:
            p("--- issuing reverse-continue #%d (this may take a while) ---" % self.step)
            gdb.execute("reverse-continue")
        except Exception as e:
            p("reverse-continue issue failed: %s" % e)
            self.finish()

    def on_stop(self, event):
        if not getattr(self, "armed", False):
            return
        p("[stop event received for step %d]" % self.step)
        try:
            gdb.execute("bt 15")
        except Exception as e:
            p("bt failed: %s" % e)
        try:
            gdb.execute("delete display")
            gdb.execute("delete")  # clears all breakpoints/watchpoints; safe here, nothing else active
        except Exception as e:
            p("cleanup failed: %s" % e)
        self.armed = False
        self.next_watch()

    def finish(self):
        p("=== REVERSE WALK FINISHED ===")
        try:
            gdb.execute("quit")
        except Exception as e:
            p("quit failed: %s" % e)

walker = ReverseWalker(30)
gdb.events.stop.connect(walker.on_stop)

gdb.execute("set pagination off")
gdb.execute("set confirm off")
gdb.execute("set breakpoint pending on")

CrashBreak("predModule.cpp:297", internal=False)
p("READY-RR-SETUP2")
