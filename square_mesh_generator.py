from sys import stdout
def generate_square_mesh(n, grain=1, f=stdout):
	'''
	n is the length of one side
	'''
	f.write('OFF\n')

	nVertices = pow(n+1, 2)
	nFaces = 2 * pow(n, 2)

	f.write('%i %i 0\n' %(nVertices, nFaces))

	for r in range(n+1):
		for c in range(n+1):
			f.write('%.2f %.2f %.2f\n' %(r * grain, c*grain, 0))

	for r in range(n):
		for c in range(n):
			bot_left = (n+1) * r + c
			bot_right = (n+1) * r + (c+1)
			top_left = (n+1) * (r+1) + c
			top_right = (n+1) * (r+1) + (c+1)
			f.write('3 %i %i %i\n' %(bot_left, top_left, top_right))
			f.write('3 %i %i %i\n' %(bot_left, bot_right, top_right))

if __name__ == '__main__':
	n = int( raw_input('How many sides would you like in your square mesh?\n') )
	generate_square_mesh(n, grain=0.01, f=open('Meshes/square/square.off', 'w'))

