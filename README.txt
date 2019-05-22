
Welcome to PVTfree
------------------

What is PVTfree?  It is an Equation of State (EoS) PVT program, that is free!  

It is intended for use by Reservoir & Petroleum Engineers for the analysis of
PVT reports and the preparation of data suitable for use in reservoir and
production simulation software.

There are NO warranties given for its use.  And support and guidance will only
be offered if agreed in advance with the original developer, Steve Furnival,
of HoBoil Limited, UK: email dev.pvtfree@gmail.com.

It is written in Python and uses modules from the NumPy, SciPy and MatPlotLib
libraries; it is expected that any potential user has enough understanding
(or who has an IT department with sufficient understanding) to 'make it work'.

At the time of writing, there is no Graphical User Interface (GUI) for PVTfree.
You must create an ASCII file of keywords and data (a bit like creating an
input 'deck' for a simulator) which you then run in the Command Prompt.
Anyone interested in developing a GUI for PVTfree?  Drop me a mail to the
address above.  I can sketch some possible windows for you.

After a training course using the code recently, I did appreciate the need
for 'help' in building the initial dataset and so wrote an Excel/VBA sheet
to aid the novice user.  Frankly, this might cause more problems than it
solves, but we'll see how it goes.

Anyone wanting to develop additional calculation modules - great, welcome on
board!  The current code is missing miscibility-type calculations (including
a 1D simulator to model slimtubes), third-phase calculations (waxes, hydrates,
asphaltenes) and much more.  Some of the existing algorithms can defintely be
improved.  I tried to develop a BFGS-based stability test but couldn't get
it to work reliably so I dropped back to the 'standard' SS/GDEM unlike the
2-phase flash where I use BFGS to finish 'hard' problems.  Also, the 'fast' 
Michelsen scheme for calculating the phase envelope is probably only half 
finished; I couldn't get the perioidic update of the saturated lines working,
nor the internal vapour fraction lines and don't ask about the critical
point calculations - they were a glorious failure in my hands.  Maybe you
can do better?

I'm reasonably confident that if you present the right data in the right way
it will work.  Problems, and I'm sure there will be many, will most likely
occur if you make a mistake in the data entry.  As a first port of call, 
look at the datasets I've provided with the release and try to make your
dataset 'look' like mine.  If you still have problems, then its time to test
your Python skills.  "But, I can't code in Python", I hear you say.  Nor could
I twelve months ago - I taught myself as I'd expect any engineer worth their
salt to do.  Remember the root of the word engineer is 'gene'.  Yes, we are
all genii, but don't tell the geoscientists, they'll get upset.
