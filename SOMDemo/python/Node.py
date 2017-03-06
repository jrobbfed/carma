
class Node(object):
	def __init__(self, numweights, x, y, left, right, top, bottom):
		self.weights = np.random.random(numweights)
		self.x = left + (right - left) / 2. 
		self.y = top + (bottom - top) / 2.