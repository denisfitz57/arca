import random, copy

INFINITY = 100000
EPSILON = .001

"""
================================================================
		 	Choose a dictionary key weighted by value
================================================================
"""
	
def weighted_choice(d):
	choices = d.items()
	total = sum(w for c, w in choices)
	if not total:
		return 0
	r = random.uniform(0, total)
	upto = 0
	for c, w in choices:
		if upto + w >= r:
			return c
		upto += w
	assert False, "Shouldn't get here"

"""
================================================================
		 	Parse strings into letters or notes
================================================================
"""
		
letterPCs = {'C': 0, 'D': 2, 'E':4, 'F':5, 'G':7, 'A':9, 'B':11}
PCletters = {0: 'C', 1: 'C#', 2: 'D', 3: 'Eb', 4: 'E', 5: 'F', 6: 'F#', 7: 'G', 8: 'G#', 9: 'A', 10: 'Bb', 11: 'B'}
def parse_lettername(s, lastOctaveSymbol = -1):
	s = s.upper()
	letterName = None
	for i, c in enumerate(s):
		letterName = letterPCs.get(c, None)
		if letterName != None:
			break
	if letterName is None:
		return parse_number(s)
	newS = s[i+1:]
	octaveSymbols = ''.join([x for x in newS if x.isdigit()])
	if not octaveSymbols:
		octave = lastOctaveSymbol
	else:
		octave = int(octaveSymbols)
	modifier = -1*(newS.count('B') + newS.count('-')) + (newS.count("#"))
	return letterName + modifier + (octave + 1)*12

def parse_number(s):
	if s.count('.') > 0:
		firstDot = s.index('.')
		newS = s[:firstDot] + '.' + s[firstDot+1:].replace('.','')
		return float(''.join([x for x in newS if x.isdigit() or x in '.-'])) 
	elif s.count('/') > 0:
		firstSlash = s.index('/')
		return float(s[:firstSlash])/float(s[firstSlash+1:].replace('/',''))
	else:
		digitChars = ''.join([x for x in s if x.isdigit() or x == '-'])
		if digitChars:
			return int(digitChars)
		return None

def spell_MIDI(mNotes = [], octaves = None):
	if type(mNotes) is int or type(mNotes) is tuple:
		mNotes = [mNotes]
	if octaves is None and max(mNotes) > 12:
		octaves = True
	return [PCletters.get(x%12, '?') + (str(-1 + int(x/12)) if octaves else '') for x in mNotes]

"""
================================================================
			NR inversions (by PCs and numerical)
================================================================
"""	
		
def NR_inversion(pitches, iList, modulus = 12):			# iList is a pair of PCs to invert around
	sortedPCs = sorted([x % modulus for x in pitches])
	tDict = NR_inversion_dict(sortedPCs, iList, modulus)
	return [x + tDict[x%modulus] for x in pitches]

def NR_inversion_by_number(pitches, inversionNumber = 0, extraRotations = 0, dualistLabels = False, modulus = 12):
	sortedPCs = sorted([x % modulus for x in pitches])
	tDict = NR_inversion_by_number_dict(sortedPCs, inversionNumber, extraRotations, dualistLabels, modulus)
	#print(pitches, tDict)
	return [x + tDict[x%modulus] for x in pitches]

def NR_inversion_dict(sortedPCs, iList, modulus = 12):
	tDict = {}
	indexNumber = sum(iList)
	targetPC = iList[0]
	pInv = sorted([(indexNumber - x) % modulus for x in sortedPCs])
	pSort = rotate_to_start(sortedPCs[:], targetPC, modulus)
	pInv = rotate_to_start(pInv, targetPC, modulus)
	if (not pSort) or (not pInv): return False
	for i, pc in enumerate(pSort):
		pc = pSort[i]
		path = pInv[i] - pc
		tDict[pc % modulus] = path
	return tDict

def NR_inversion_by_number_dict(sortedPCs, inversionNumber = 0, extraRotations = 0, dualistLabels = False, modulus = 12):
	""" inversion by number, either dualist or nondualist """
	
	gnf = geometrical_normal_form_local(sortedPCs, invert = dualistLabels, extraRotations = extraRotations, modulus = modulus)
	cardinality = len(gnf)
	
	if cardinality % 2 == 0:
		v1 = int(-inversionNumber / 2) % cardinality
		v2 = (v1 + 1) % cardinality
		if inversionNumber % 2 == 1:
			v1 = (v1 - 1) % cardinality	
	else:
		v1 = int((-inversionNumber*(len(gnf)+1)/2) % cardinality)
		v2 = (v1 + 1) % cardinality
	
	transp = None
	
	for i in range(modulus):
		if all([(x + i) % modulus in sortedPCs for x in gnf]):
			transp = i
			break
			
	if transp is not None:
		note1 = (gnf[v1] + transp) % modulus
		note2 = (gnf[v2] + transp) % modulus	
	else:
		if not dualistLabels:
			return False	
		index = None
		for i in range(modulus):
			if all([(i - x) % modulus in sortedPCs for x in gnf]):
				index = i
				break
		if index is None:
			return False
		note1 = (index - gnf[v1]) % modulus
		note2 = (index - gnf[v2]) % modulus
		
	return NR_inversion_dict(sortedPCs, [note1, note2], modulus = modulus)

def rotate_to_start(l, targetNote = 0, modulus = 12):
	if l[0] == targetNote:
		return l
	newList = l
	for i in range(len(l)):
		newList = newList[1:] + [newList[0] + modulus]
		if newList[0] == targetNote:
			return newList
	return False

"""
================================================================
				Put a set into normal form
================================================================
"""
		
def geometrical_normal_form_local(inList, invert = False, modulus = 12):				
	
	"""
	
	geometrical normal form is different from standard normal form insofar as it puts the smallest interval in the first position
	
	"""
	
	theLen = len(inList)
	listOfPCs = sorted([(k - inList[0]) % modulus for k in inList])
	listOfInts = [listOfPCs[i] - listOfPCs[i-1] for i in range(1, theLen)]
	listOfInts.append(modulus - listOfPCs[-1])
	currentBest = listOfInts[:]
	
	newChallenger = currentBest
	for i in range(1, theLen):
		newChallenger = newChallenger[1:] + [newChallenger[0]]
		if newChallenger < currentBest:
			currentBest = newChallenger
			
	if invert:	
		newChallenger = [currentBest[0]] + currentBest[-1:0:-1]
		if newChallenger < currentBest:
			currentBest = newChallenger
			
		for i in range(1, theLen):
			newChallenger = newChallenger[1:] + [newChallenger[0]]
			if newChallenger < currentBest:
				currentBest = newChallenger
	
	if extraRotations:
		minInt = currentBest[0]
		numInts = currentBest.count(minInt)
		extraRotations = extraRotations % numInts
		i = 0
		count = 0
		while count < extraRotations:
			i = (i + 1) % len(currentBest)
			if currentBest[i] == minInt:
				count += 1
		currentBest = currentBest[i:] + currentBest[:i]
	
	outList = [0]
	for i in currentBest[:-1]:
		outList.append(outList[-1] + i)
	return outList

def get_chord_element_names(PCs, gnf, modulus = 12):
	"""	identifies chord elements by their normal forms; 
		returns a dictionary mapping the indices of PCs to those of gnf, 
		plus True or False depending on whether there is inversion involveds"""
	chordNames = {}

	transp = None
	
	inverted = False

	for i in range(modulus):
		working = [(x + i) % modulus for x in gnf]
		if all([x in PCs for x in working]):
			transp = i
			break

	if transp is None:
		index = None
		for i in range(modulus):
			working = [(i - x) % modulus for x in gnf]
			if all([x in PCs for x in working]):
				index = i
				break
		inverted = True
		if index is None:
			return {}, False
				
	for i, pc in enumerate(PCs):
		target = working.index(pc)
		chordNames[i] = target
		working[target] = None

	return chordNames, inverted


"""
================================================================
			Label a voiceleading by transformations
================================================================
"""
	
class VLLabel():
	
	"""
	
	Label voiceleadings: input two chords of pitches and it will return a string of voice leading transformations generating the VL 
	
		Tx = chromatic transposition
		tx = transposition along the intrinsic scale
		cx = voice crossing of notes -x and -x + 1 in the normal form
		Jx = neo-Riemannian voice leading (dualist label)
		Kx = neo-Riemannian voice leading (nondualist label)
		P(x, y, z, ...) = perturbation moving voice 1 by x, voice 2 by y, voice 3 by z ...
	
	When findInversions is true the labels can include NR voice leadings among its transpositions,
		otherwise inversionally related chords are treated as distinct
	
	When dualistLabels is True, self.invert is true by default
		if dualistLabels is True, and startChord is not in the inversional normal form region, then the whole voice leading is inverted
	
	When factorCrossings is True crossings Cx are factored as txC0t-x.
	
	TODO: 
		needs testing: make sure everything is consistent with the labels used by arca.transform
	
	"""
	
	def __init__(self, startPos = [60, 64, 67], endPos = [67, 64, 60], findInversions = True, dualistLabels = True, factorCrossings = False, autoPrint = True, modulus = 12):
		
		self.originalLastNotes = startPos
		self.originalCurrentNotes = endPos
		
		self.lastNotes = startPos
		self.currentNotes = endPos
		
		self.modulus = modulus
		self.chordLen = len(self.lastNotes)
		
		self.transpositionRegion = factorCrossings
		self.invert = findInversions
		self.dualistLabels = dualistLabels
		
		if dualistLabels:
			self.invert = True
		
		self.iLabel = None
		self.newSet = False
		self.gnf1invert = None
		self.gnf2invert = None
		self.normal_form_VL = False
		self.intrinsicTransp = 0
		self.totalTransp = 0
		self.dualistFlag = False
		
		self.VLstring = self.identify_voiceleading()
		
		if autoPrint:
			print(self.VLstring)
	
	def make_labels(self):
		
		self.gnf1noinvert = geometrical_normal_form_local(self.lastNotes, invert = False, modulus = self.modulus)
		
		"""if we are using dualist labels, we invert the whole voice leading if the original chord is not in inversional normal form"""
		if self.dualistLabels:
			self.gnf1invert = geometrical_normal_form_local(self.lastNotes, invert = True, modulus = self.modulus)
			if self.gnf1invert != self.gnf1noinvert:
				self.lastNotes = [120 - x for x in self.lastNotes]
				self.currentNotes = [120 - x for x in self.currentNotes]
				self.gnf1noinvert = self.gnf1invert
				self.dualistFlag = True									# Flag to note we have inverted the voice leading
		
		self.set1 = [x % self.modulus for x in self.lastNotes]
		self.chordNames1, junk = get_chord_element_names(self.set1, self.gnf1noinvert, self.modulus)
		
		self.gnf2noinvert = geometrical_normal_form_local(self.currentNotes, invert = False, modulus = self.modulus)
		
	def identify_voiceleading(self):
		
		self.make_labels()
		self.get_crossing_times()
		self.crossingLabels = self.calculate_crossings()
		
		"""case 1: they are related by transposition"""
		if self.gnf1noinvert == self.gnf2noinvert:
			self.totalTransp, self.intrinsicTransp = self.find_intrinsic_transposition(self.lastNotes, self.currentNotes)
			return self.formulate_VL_string()
		elif self.dualistLabels or self.invert:
			"""case 2: they are related by inversion"""	
			if not self.gnf1invert:
				self.gnf1invert = geometrical_normal_form_local(self.lastNotes, invert = True, modulus = self.modulus)
			self.gnf2invert = geometrical_normal_form_local(self.currentNotes, invert = True, modulus = self.modulus)
			if self.gnf1invert == self.gnf2invert:
				self.totalTransp, self.iLabel = self.find_NR_inversion(self.lastNotes, self.currentNotes)
				return self.formulate_VL_string()
		
		"""case 3: unrelated sets; we have to find a perturbation that turns normal form into the other."""
		
		"""
		If self.dualistLabels is True than chord1 is in the inverted normal form region, so chord2 should be as well;
		if it is False then we are in the transposition region, and so we do not want to invert at all.
		
		"""
						
		if self.dualistLabels:
			self.gnf2 = self.gnf2invert
		else:
			self.gnf2 = self.gnf2noinvert
		
		#TODO: this code seems no longer to be necessary
		
		"""startChord is not in the inversional normal-form region but endChord is"""
		"""
		if self.dualistLabels and (self.gnf1invert != self.gnf1noinvert) and (self.gnf2invert == self.gnf2noinvert):
			self.gnf2 = geometrical_normal_form_local([-x % 12 for x in self.currentNotes], invert = False, modulus = self.modulus)"""
		
		"""calculate the perturbation, and work with the perturbed form"""
		self.normal_form_VL = [self.gnf2[i] - self.gnf1noinvert[i] for i in range(self.chordLen)]
		self.adjustedChord1 = [self.lastNotes[i] + self.normal_form_VL[self.chordNames1[i]] for i in range(self.chordLen)]
		
		"""case 3a: perturbed start chord and end chord are related by intrinsic transposition"""
		self.totalTransp, self.intrinsicTransp = self.find_intrinsic_transposition(self.adjustedChord1, self.currentNotes)
		
		"""case 3b: perturbed start chord and end chord are related by neo-Riemannian inversion"""
		if self.intrinsicTransp is None:
			self.totalTransp, self.iLabel = self.find_NR_inversion(self.adjustedChord1, self.currentNotes)
		
		return self.formulate_VL_string()
	
	def formulate_VL_string(self):
		vlString = ''
		if self.normal_form_VL:
			vlString = f"P{[self.normal_form_VL[self.chordNames1[i]] for i in range(self.chordLen)]}".replace('[', '(').replace(']', ')').replace(' ', '')
		if self.crossingLabels:
			if self.transpositionRegion:
				cString = [f't{x}C0' for x in self.transposition_region_crossings()]
				cString = ''.join(cString)[:-2]								# strip off last c0
				vlString += cString
			else:
				vlString += ''.join(['C'+str(x) for x in self.crossingLabels])				
		if self.intrinsicTransp:
			vlString += f't{self.intrinsicTransp}'
		if self.iLabel is not None:
			if self.dualistLabels:
				vlString += f'J{self.iLabel}'
			else:
				vlString += f'K{self.iLabel}'
		if self.totalTransp:
			vlString += f'T{self.totalTransp}'.replace('.0', '')
		return vlString
		
	def transposition_region_crossings(self):
		tList = []
		lastTransposition = 0
		for c in self.crossingLabels:
			tList.append((lastTransposition + c)%self.chordLen)
			lastTransposition = -c
		tList.append((lastTransposition + self.intrinsicTransp) % self.chordLen)
		return tList
	
	def find_NR_inversion(self, startNotes, endNotes):
		
		endSum = sum(endNotes)
		endPCs = sorted([round(x % self.modulus, 2) for x in endNotes])
		
		foundChord = False
		
		for i in range(self.chordLen):
			newChord = NR_inversion_by_number(startNotes, i, dualistLabels = self.dualistLabels, modulus = self.modulus)
			tLevel = (endSum - sum(newChord))/self.chordLen
			transposedNewPCs = [round((x + tLevel)%self.modulus, 2) for x in newChord]
			if all([x in endPCs for x in transposedNewPCs]):
				foundChord = True
				break
		
		if foundChord:
			return tLevel, i
			
		return 0, None
		
	def find_intrinsic_transposition(self, startChord, endChord):
		
		sumDiff = sum(endChord) - sum(startChord)
		avgSumDiff = sumDiff/self.chordLen
	
		transposedChord1 = [startChord[i] + avgSumDiff for i in range(self.chordLen)]
		
		endPCs = sorted([round(x % self.modulus, 2) for x in endChord])
		startPCs = sorted([round(x % self.modulus, 2) for x in transposedChord1])
		avgSum = self.modulus/self.chordLen
		
		for i in range(self.chordLen):
			if all([startPCs[i] == endPCs[i] for i in range(self.chordLen)]):
				return -(-sumDiff + (i * self.modulus))/self.chordLen, i
			startPCs = sorted([round((x - avgSum)%self.modulus, 2) for x in startPCs])
		
		return 0, None	
		
	def calculate_crossings(self):
		
		"""NB: chordNames is a local variable since we relabel as we go"""
	
		chordNames = {x[0]:x[1] for x in self.chordNames1.items()}
		lastCrossing = None
		crossingLabels = []
	
		for t in sorted(self.crossingTimes.keys()):
			for v1, v2, distance in self.crossingTimes[t]:
				voice1 = chordNames[v1]
				voice2 = chordNames[v2]
			
				if self.chordLen != 2:
					if voice2 == (voice1 + 1) % self.chordLen:
						lastCrossing = (-voice1)%self.chordLen
					elif voice1 == (voice2 + 1)%self.chordLen:
						lastCrossing = (-voice2)%self.chordLen
					else:
						print('Problem labeling crossings:', voice1, voice2)
						continue
				else:
					if lastCrossing == None:
						if distance > self.modulus/2:
							lastCrossing = 1
						else:
							lastCrossing = 0
					else:
						lastCrossing = (lastCrossing + 1)%2
				crossingLabels.append(lastCrossing)
					
				chordNames[v1] = voice2
				chordNames[v2] = voice1
			
		return crossingLabels
	
		
	def get_crossing_times(self):
		self.moveVector = [self.currentNotes[i] - self.lastNotes[i] for i in range(self.chordLen)]
		self.crossingTimes = {}
	
		for firstVoice in range(self.chordLen - 1):
			for secondVoice in range(firstVoice + 1, self.chordLen):
				relativeVelocity = self.moveVector[secondVoice] - self.moveVector[firstVoice]
				distance = (self.lastNotes[secondVoice] - self.lastNotes[firstVoice]) % self.modulus
				if relativeVelocity > 0:
					distance = self.modulus - distance
				elif relativeVelocity == 0:
					continue
				vel = abs(relativeVelocity)
				crossingTime = distance/vel
				while crossingTime < 1:
					self.crossingTimes.setdefault(crossingTime, []).append((firstVoice, secondVoice, distance))
					crossingTime += self.modulus/vel
	
	def __str__(self):
		return self.VLstring
		
"""
================================================================
			Find the minimum VL between two sets
================================================================
"""
		
def minimum_vl(first, second, sort = True, modulus = 12):		# do better for chords with different lengths
	minimum_vl.fullList = []
	if len(second) > len(first):
		newNotes = []
		for i in range(0, len(second) - len(first)):
			newNotes.append(first[random.randrange(0, len(first))])
		first = first + newNotes
	elif len(second) < len(first):
		for i in range(0, len(first) - len(second)):
			first.pop(random.randrange(0, len(first)))
	firstPCs = sorted([p % modulus for p in first])
	secondPCs = sorted([p % modulus for p in second])
	secondPCs = secondPCs[1:] + [secondPCs[0] + modulus]
	currentBest = []
	currentBestSize = INFINITY									# very large number
	for i in range(0, len(firstPCs) + 1):
		newSize = simple_vl_size(firstPCs, secondPCs, modulus = modulus)
		newPaths = [[firstPCs[i], secondPCs[i] - firstPCs[i]] for i in range(len(firstPCs))]
		minimum_vl.fullList.append([newPaths, newSize])
		if newSize < currentBestSize:
			currentBestSize = newSize
			currentBest = newPaths
		secondPCs = [secondPCs[-1] - modulus] + secondPCs[:-1]
	minimum_vl.size = currentBestSize
	if sort:
		minimum_vl.fullList = sorted(minimum_vl.fullList, key = lambda x: x[1])
	return currentBest
	
"""
================================================================
				Flatten a list containing sublists
================================================================
"""

"""
	the ignore function describes lists that should be preserved	
	
"""
	
def default_ignore_function(item):
	if type(item[0]) is str:
		return True
	return False
	
def flatten_note_list(l, ignoreFunction = default_ignore_function):
	
	flatList = []
	noteGroups = []
	
	for item in l:
		curLen = len(flatList)
		if (type(item) is not list) or ignoreFunction(item):
			noteGroups.append(curLen)
			flatList.append(item)
		else:
			item = flatten_list(item, ignoreFunction = ignoreFunction)
			noteGroups.append(slice(curLen, curLen + len(item)))
			flatList += item
			
	return flatList, noteGroups

def flatten_list(l, ignoreFunction = default_ignore_function):
	out = []
	for item in l:
		if (type(item) is not list) or ignoreFunction(item):
			out.append(item)
		else:
			out += flatten_list(item)
	return out

"""
================================================================
		Group a split-up list of strings by parentheses
================================================================

   Turns ['ca', 'p(', 'y', ')', 'bara'] into ['ca', p(y)', 'bara']
	
"""


def parse_parentheses(sSplit, joinSymbol = ' ', pairs = ['{}', '()']):
	for openItem, closeItem in pairs:
		sSplit = parse_paren_unit(sSplit, openItem, closeItem, joinSymbol)
	return sSplit

def parse_paren_unit(sSplit, openItem = '{', closeItem = '}', joinSymbol = ' '):
	out = []
	targetList = out
	depth = 0
	
	for i, s in enumerate(sSplit):	
		delta = s.count(openItem) - s.count(closeItem)
		if depth == 0 and delta > 0:
			targetList = []
			
		targetList.append(s)
		
		if depth > 0 and delta == -depth:
			out.append(joinSymbol.join(targetList))
			targetList = out
			
		depth = depth + delta
	
	return out
