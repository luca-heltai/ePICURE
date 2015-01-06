def bold(msg):
	return u'\033[1m%s\033[0m' % msg

def color(this_color, string):
	"""There are 8 colors, ANSI codes 30 to 37, which can have the bold modifier, making 16 colors.
	>>> for i in range(30, 38):
	>>> c = str(i)
	>>> print('This is %s' % color(c, 'color ' + c))
	>>> c = '1;' + str(i)
	>>> print('This is %s' % color(c, 'color ' + c))
	"""
	return "\033[" + this_color + "m" + string + "\033[0m"