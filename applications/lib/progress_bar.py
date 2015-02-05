from lib.font import *
import sys
import fcntl
import termios
import struct

class progress_bar(object):

	def __init__(self, tot=100, lenght=10):
		self.cp='/-\|'
		self.bar_lenght = lenght
		self.tot = tot

	def startprogress(self, title):
		"""Creates a progress bar 40 chars long on the console
		and moves cursor back to beginning with BS character"""
		sys.stdout.write(title + ": [" + "-" * self.bar_lenght + "]" + chr(8) * (self.bar_lenght+1))
		sys.stdout.flush()

	def progress(self, x):
		"""Sets progress bar to a certain percentage x.
		Progress is given as whole percentage, i.e. 50% done
		is given by x = 50"""
		y = int(x)%4
		z = int((x/float(self.tot))*self.bar_lenght)                      
		sys.stdout.write("#" * z + self.cp[y] +"-" * (self.bar_lenght-1 - z) + "]  "+ bold(str(int(x))+"/"+str(self.tot)) + chr(8) * (self.bar_lenght+4+len(str(int(x)))+len(str(self.tot)) ))
		sys.stdout.flush()

	def endprogress(self):
		"""End of progress bar;
		Write full bar, then move to next line"""
		sys.stdout.write("#" * self.bar_lenght + "]\n")
		sys.stdout.flush()

class all_line_progress_bar(object):

	def __init__(self):
		self.COLS = struct.unpack('hh',  fcntl.ioctl(sys.stdout, termios.TIOCGWINSZ, '1234'))[1]

	def progress(self,current, total):
		prefix = '%d / %d' % (current, total)
		bar_start = ' ['
		bar_end = '] '
		bar_size = self.COLS - len(prefix + bar_start + bar_end)
		amount = int(current / (total / float(bar_size)))
		remain = bar_size - amount
		bar = '#' * amount + ' ' * remain
		return bold(prefix) + bar_start + bar + bar_end

	def bar(self, current, total):
		sys.stdout.write(self.progress(current,total) + '\r')
		sys.stdout.flush()