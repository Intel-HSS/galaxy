import sys

file_to_read = sys.argv[1]
file_to_write = sys.argv[2]

mapped_chars = { '__gt__': '>' ,
                 '__lt__': '<',
                 '__sq__': "'" ,
                 '__dq__': '"',
                 '__ob__': '[',
                 '__cb__': ']',
                 '__oc__': '{',
                 '__cc__': '}',
                 '__at__': '@',
                 '__cn__': '\n',
                 '__cr__': '\r',
                 '__tc__': '\t',
                 '__pd__': '#'}

def read_and_write_file(read_this, write_this):
	f1 = open(read_this, 'r')
	lines = f1.read()
	for i, j in mapped_chars.iteritems():
		lines = lines.replace(i, j)
	f2 = open(file_to_write, 'w')
	f2.write(lines)
	f2.close

read_and_write_file(file_to_read, file_to_write)
