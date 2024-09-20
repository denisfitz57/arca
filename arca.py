import random, copy
from arcatools import *				# contains utility routines used throughout
import DTXML

DEFAULTDIRECTORY = '/Users/dmitri/Source Code/Python/Arca Programs/'
DEFAULTOUTPUT = '/Users/dmitri/Desktop/arcaoutput.musicxml'

"""
Arca, a musical programming language (after Athanasius Kircher's Arca Musicarithmica, described in his 1650 book Musica Universalis) 

VERSION 9/20/2024

	Fundamental classes:

	User:
		Scale					Nestable scales
		ScalePattern			A concrete collection of notes [replaces the more awkward Chord/Arpeggiation pair]
		
		To be phased out:
			Chord					A concrete collection of notes, can have a default Arpeggiation, contains a Scale object for calculations
			Arpeggiation			A "structured arpeggiation" that groups and orders a chord's notes, including nonharmonic tones
		
		Pattern					Combines arpeggiations from multiple chords into a flexible compositional object (TODO)
		Program					Text interface to Arca
		VelocityGenerator		simple routines for creating velocity structures (for MIDI applications)
		RhythmGenerator			simple routines for creating rhythmic patterns [TODO: add recursive nestability]
	
	Internal:
		CoreVoiceLeading		Voice-leading routines common to Chord and Scale, not for direct use [but start here I think]
									kind of pointless now that I have decided to phase out Chord, could be merged into Scale

		Timeline				Used by Program, executes Arca functions in sequence
		RepeatFunctionWrapper	Silly little class to allow indexing of repeatFunctions
	
	Implement lookforward for forward neighbors?
		in NHT, have certain notes wait until later notes are calcuated

	Motive copying
	
	Add a FILENAME attibute that determines the filename on save?
	
	SMALL AND IMPORTANT:
		add a more general undo transformation
	
	MEDIUM:
		Make a superclass of Arpeggiation that 
	 		1. supports VelocityGenerators
			2. supports beat-based rhythms (that is, you call it every beat and sometimes get notes)
			3. has the ability to automatically apply the chord's next() method after output
	
	TODO IMPORTANT:
		In brython False == None is True, whereas None == False is False (I believe javascript converts to the first type)?  

		Change spelling elegantly and quickly, at anypoint in the stream (sharps, flats)
			good spelling is important and hard

		griegNOWEB failed because of the @m syntax, but only in brython???  That is weird.  It works fine when converted to beats.
			- solved by converting @m to integers, when they should be converted.  Should this happen elsewhere?
			- and why was this making problems?  Makes me think the beat regularization is error-filled.  Check that out more thoroughly?

		Chrome won't play Spiegel.  Too much audio info cued too far in advance?

	TODO:
		Rethink chord -- a chord should be a subclass of scale, or should include a scale in it
		NHT syntax: 
			level 0 current chord
			level 1-n containers
			level None = chromatic
	
	TODO:

		1. web interface:
			test microtonality

		2. spelling

		3. Build up Pattern's capacities
				pattern should be able to layer two objects on top of each other

		4. A Process object?  Or can chord do it?

		5. Test and expand Program

		6. Incorporate PitchMotive's capacities into Arpeggiation and Pattern;
				maybe a superclass of Apreggiation and Pattern supporting VelocityGenerators, RhythmGenerators, and groups

		7. Build a test suite

	Smaller things:

		1. allow manual control and initialization of multisets in Chord, using four digits

		2. generalize repeat_function_permutation, so that repeat_functions are only a special case
			I think find_transposition_and_permutation should do this, it is just a matter of consolidating

		3. eliminate redundant or outdated methods, streamline

		4. Incorporate range checking into objects?

		Replace the arpeggiation "chain" keyword with a more general "linkinterval" that connects the last note of the previous pattern to the first note of the second.

		Allow timelines to play notes from a scale?  Can this be done without scale degrees?  Should I use variables like Scalename[X = C4], so then later you can write Scalename[x]?

		Allow arpeggiation patterns to be instantiated with chordal intervals like (+1, +2) -- maybe this can be done currently with chordal skips?

		Test @m, for referring to time by measure numbers

	Philosophical principles:
		- use pure python
		- avoid using scale degree numbers when possible (everything is initialized with MIDI)
		- constructing objects is often easier with text
		- use python for anything complicated; use the text language for simple things

"""

class Default():
	
	def set_default_arguments(self, kwargs):
		self.defaultArgs = copy.deepcopy(type(self).defaults)
		for k, val in kwargs.items():
			self.defaultArgs[k] = val
		
		for k, val in self.defaultArgs.items():
			setattr(self, k, val)
		
class CoreVoiceleading(Default):
	
	"""
	
	Voice leading routines shared by both Scale and Chord; the philosophy is to put as much here as possible

		get_transpositions:	find the basic voice leading for a given chord and scale
		double_transpose:	simultaneous transposition along input and output degrees
	
	Voice leading changes should go through appply_VL_dict
	
	Originally the thought was that a Chord would not have a scale object attached; I later decided a default Scale provided some benefits

	We assume:
		the object has a currentNotes attribute
			- for a chord, these are just its 'raw' or harmonic notes, shorn of duplications
			- for a scale, a 'reference octave' of scale degrees
		the object has a lastNotes attribute
		the object has PCs (should be in the same order as currentNotes)
		the object has sortedPCs	(sorted version of PCs)
	
		the object has an .update method that is called upon small changes to its notes (i.e. non-cardinality-changing)
		the object may or may not have a container, which is a scale
		the object's notes can be indexed with getitem
	
		output_note is used for multisets; 
			multisets have duplicate scale degrees represented by the fourth decimal place:
				60.0000			-> first scale degree = 60
				60.0001			-> second scale degree = 60
			output_note quanties to cents, sending both of these to 60.00
	
		repeat_functions are an easy way of creating chord progressions; the .next method cycles through these
	
	Main methods:
		transform -- apply voice leading transformations via text
		lattice -- another set of voice leading transformations based on scrambling the circle of fifths
	
	"""					
	
	transformSubstitutions = {'VOiCELEAD': 'M', 'VL(': 'M(', 'VOiCEEXCHANGE': 'C', 'VE(': 'C(',}
	
	def initialize_core(self):
		"""create datastructures that are assumed to exist"""
		self.permDict = {}
		self.lastVoiceleading = {}
		self.lastNotes = []
		self.suppressChangeFunc = False
		self.changeFunc = False
		self.lastTransposition = 0
		self.repeatFunctionWrapper = RepeatFunctionWrapper(self, [])
		self.progression = []
		self.iterations = -1
		
	def get_basic_VL(self):
		vl = CoreVoiceleading.get_transpositions(len(self.currentNotes), self.modulus)			# got these reversed, VL1 should be the sum-changing one
		
		if len(vl) > 1:
			self.basicVL = vl[1]
			self.basicVL2 = vl[0]
		else:
			self.basicVL = vl[0]
			self.basicVL2 = []
		self.basicVLSumChange = self.basicVL[0]*len(self.currentNotes) + self.basicVL[1] * self.modulus
	
	def apply_VL_dict(self, vlDict, globalTransposition = 0):
		self.lastVoiceleading = vlDict
		self.lastNotes = self.currentNotes
		self.currentNotes = [x + vlDict.get(x%self.modulus, 0) + globalTransposition for x in self.currentNotes]
		self.update()
		if self.changeFunc and (not self.suppressChangeFunc):
			self.changeFunc(self)
	
	def output_MIDI_note(self, n):
		if self.container:
			return self.output_note(self.container.__getitem__(n, recurse = True)) 
		return self.output_note(n)
		
	def output_note_regular(self, n):
		return n
	
	def output_note_multiset(self, n):
		if type(n) is int:
			return n
		n = round(n, 2)
		if n == int(n):
			return int(n)
		return n
	
	def raw_content(self):
		return [self.chromaticPC(x) for x in self.currentNotes]	
	
	def chromaticPC(self, scaleDegree):
		return (self.container.__getitem__(scaleDegree, recurse = True) if self.container else scaleDegree) % 12
			
	def Tt(self, *args):
		self.double_transpose(*args)
	
	def voicelead(self, thePairs):
		"""
		pass a list of [pitch, path] pairs for a voice leading at the scalar level
		Warning: NO ERROR CHECKING! if your PC isn't in the scale, nothing happens
		
		usually, you will do this at the 60+x pitch level: for 0, 2, 4 in the diatonic (C, E, G) use (60, 62, 64)
		"""
		self.apply_VL_dict({x[0]%self.modulus:x[1] for x in thePairs})
	
	def apply_basicVL(self, applications = 0):
		self.double_transpose(self.basicVL[0]*applications, self.basicVL[1]*applications)
	
	def multischema(self, listOfVLLists):
		"""
		path a schematic voice leading like [[0, -1], [4, 0], [7, 1]] for (C, E, G)->(B, E, G#)
		this routine will apply this voice leading schema to the current chord, assuming it is a transposition of the input chord
		"""
		for VLList in listOfVLLists:
			result = self.transpositional_VL_schema(VLList)
			if result: return True
		return False
	
	def reverse_schema(self, schema):
		if type(schema) is dict:
			schema = list(schema.items())
		return [[sum(x)%self.modulus, -x[1]] for x in schema]
	
	def invert_schema(self, schema):
		if type(schema) is dict:
			schema = list(schema.items())
		return [[-x[0]%self.modulus, -x[1]] for x in schema]
	
	def transpositional_VL_schema(self, VLlist):
		"""
		path a schematic voice leading like [[0, -1], [4, 0], [7, 1]] for (C, E, G)->(B, E, G#)
		this routine will apply this voice leading schema to the current chord, assuming it is a transposition of the input chord
		"""
		if type(VLlist) is dict:
			VLlist = list(VLlist.items())
		foundTransp = False
		for i in range(self.modulus):
			tempChord = [(x[0] + i) % self.modulus for x in VLlist]
			if all([x in self.PCs for x in tempChord]):
				foundTransp = True
				break
		if foundTransp:
			self.apply_VL_dict({tempChord[i]:VLlist[i][1] for i in range(len(tempChord))})
			return True
		return False
	
	def minimum_VL(self, theChord = [], midiNotes = []):
		"""
		quick and dirty routine for calculating a small voice leading to any set of notes
		
		theChord is a destination chord described using scale degrees; midiNotes is a collection of midiNotes
		
		"""
		if midiNotes:
			if self.container:
				theChord = [self.container.scalarPC(x) for x in midiNotes]
			else:
				theChord = [x % 12 for x in midiNotes]
		# WARNING: cutting out this next batch of code because it doesn't make sense to me
		"""else:
			theChord = [(x + self.keyOffset)%self.modulus for x in theChord]"""					
		vl = minimum_vl(self.sortedPCs, theChord, sort = True, modulus = self.modulus)
		vlDict = {x[0]:x[1] for x in vl}
		self.apply_VL_dict(vlDict)
		
	def random_slide(self, interval = -1):
		"""randomly shift a scale note by interval"""
		availablePCs = [x for x in self.PCs if (x + interval) % self.modulus not in self.PCs]
		if not availablePCs:
			return False
		else:
			PC = random.choice(availablePCs)
			vlDict = {PC: interval}
		self.apply_VL_dict(vlDict)
		return True
	
	"""These next functions allow for the construction of progressions, useful for creating chord loops"""
	
	@property
	def progression(self):
		return self.repeatFunctionWrapper
	
	@progression.setter
	def progression(self, funcList):
		self.repeatCallChain = []
		
		for theFuncs in funcList:
			self.repeatCallChain.append(self.translate_repeat_function(theFuncs))
			
	def translate_repeat_function(self, theFuncs):
		callList = []
		
		if theFuncs:
			t = type(theFuncs[0]) 
			
			if t is not list and t is not tuple:
				theFuncs = [theFuncs]
				
			for f in theFuncs:
				
				if type(f[0]) is str:
					if type(f) is str:
						funcName = f
					else:
						funcName = f[0]
					if not (hasattr(self, funcName)):									# pass a text transformation
						callList += list(self.transform(funcName, getCallChain = True))
						continue
						
					func = getattr(self, funcName)
				else:
					func = f[0]
					
				args = []
				kwargs = {}
				if len(f) > 1:
					for i, item in enumerate(f):
						if i == 0: continue
						if type(item) is list:
							args = item
						elif type(item) is dict:
							kwargs = item
				callList.append([func, args, kwargs])
				
		return callList
	
	def apply_repeat_function(self):
		
		self.iterations += 1
		if not self.progression:
			return
		self.apply_call_chain(self.repeatCallChain[self.iterations % len(self.progression)])
		
	def next(self):
		self.apply_repeat_function()
	
	def undo(self):
		"""go back to the previous chord"""
		temp = self.lastNotes
		self.lastNotes = self.currentNotes
		self.currentNotes = temp
		self.update()
	
	def save_notes(self):
		self.initialNotes = self.currentNotes[:]
	
	def restore_initial_notes(self):
		self.currentNotes = self.initialNotes[:]
		self.update()
	
	def get_chord_from_string(self, s):
		MIDIpitches = [parse_lettername(x) for x in s.split()]
		if self.container:
			return [self.container.scale_degree(x) for x in MIDIpitches]
		return MIDIpitches
		
	def chords_to_transformations(self, chordList, schema = False):
		
		"""
		
		pass a list of midiNotes, chordlist = [...], and the program will automatically generate a progression from there; 
		used for quick initialization of rounds and sequences
		
		When schema is False you will get specific voice leadings and the final chord should have the same PCs as the first
		When it is true the progression will be a transposition VL schema and you can make a sequence
		
		"""
		
		if schema:
			funcName = 'transpositional_VL_schema'
		else:
			funcName = 'apply_VL_dict'
		
		chordList = [self.get_chord_from_string(x) if type(x) is str else x for x in chordList]
		
		l = len(chordList[0])
		prog = []
		
		for i, obj in enumerate(chordList):
			if i == 0:
				continue
			pathDict = {(chordList[i-1][j])%self.modulus:(chordList[i][j] - chordList[i-1][j]) for j in range(l)}
			prog.append([funcName, [pathDict]])
		
		return chordList[0], prog
		
	"""
		
	these next routines store the scale's PCs as MIDI notes, and then restore the scale's PCs from those notes
	this permits changing the modulus of the underlying scale
	
	"""
	
	def get_MIDI(self):
		self.lastMIDI = [self.container.__getitem__(x, recurse = True) for x in self.currentNotes] if self.container else self.currentNotes[:]
	
	def restore_from_MIDI(self, modulus = None, quantize = 0):
		if self.container:
			self.currentNotes = [self.container.scale_degree(x, direction = quantize, recurse = True) for x in self.lastMIDI]
			self.modulus = modulus if modulus else len(self.container)
		else:
			self.currentNotes = self.lastMIDI[:]
			self.modulus = modulus if modulus else 12
		self.update()
		self.get_basic_VL()
		
	def uncross(self):
		"""uncross the underlying scale
		COMPLETELY UNTESTED!
		"""
		self.currentNotes = sorted(self.currentNotes)
		self.update()
				
	"""========= label the voice leading connecting one chord to another ========="""
		
	def find_voiceleading(self, startNotes = None, endNotes = None, **kwargs):
		if not startNotes:
			startNotes = self.lastNotes
		if not endNotes:
			endNotes = self.currentNotes
		self.VLLabel = VLLabel(startNotes, endNotes, modulus = self.modulus, autoPrint = False, **kwargs)
		return str(self.VLLabel)
	
	def get_spacing(self, notes = None, startNote = None):
		if notes is None:
			notes = self.currentNotes
		if startNote is None:
			startNote = min(notes)
		out = []
		startNoteDegree = self.sortedPCs.index(startNote % self.modulus)
		for x in notes:
			if x >= startNote:
				out.append((self.sortedPCs.index(x % self.modulus) - startNoteDegree) % len(self.currentNotes) + int((x - startNote)/self.modulus)*len(self.currentNotes))
			else:
				out.append(-((startNoteDegree - self.sortedPCs.index(x % self.modulus)) % len(self.currentNotes) + int((startNote - x)/self.modulus)*len(self.currentNotes)))
		return out
		
	def apply_spacing(self, spacing = None, startNote = None):
		if not spacing: return False
		if not startNote:
			startNote = min(self.currentNotes)
		sd = self.scale.scale_degree(startNote, recurse = False)
		outDegrees = [sd + x for x in spacing]
		return [self.scale.__getitem__(x, recurse = False) for x in outDegrees]
	
	"""========= voice leading routines using the generalized circle of fifths ========="""		
	def lattice(self, vectorList = [1], generatingVL = None, sumCorrection = False):
		
		"""
		scrambles the voice leadings on the circle of fifths [= basicVL, or generatingVL], which are numbered [0, 1, 2, ...]
		a list like [1, 2] will apply vl1, then vl2.	
		
		this gives generalized melodic and harmonic minor collections
		
		sumCorretion finds the chord with the same sum as the current chord
			CURRENTLY DOES NOT WORK WITH GENERATINGVL, needs correction
		
		"""
		
		self.lastLatticeParameters = [vectorList, generatingVL, sumCorrection]
		
		if generatingVL is None:
			generatingVL = self.basicVL[:]
		
		if type(vectorList) is int:
			vectorList = [vectorList]
			
		if vectorList[0] < 0:
			vectorList = [-x for x in vectorList]
			generatingVL = [-x for x in generatingVL]
			
		if sumCorrection and self.basicVLSumChange:
			self.latticeMoves = -len(vectorList)
			self.double_transpose(self.basicVL[0]*self.latticeMoves, self.basicVL[1]*self.latticeMoves)
		else: 
			self.latticeMoves = 0
			
		self.latticeVL = generatingVL
		changeFuncStatus = self.suppressChangeFunc
		self.suppressChangeFunc = True
		vlPaths = []
		
		"""we will need this to restore the voices to their original state"""
		tempPCs = self.currentNotes[:]
		
		"""calculate the effect of the circle of fifths on the actual voices, getting as many as we need"""
		for i in range(max(vectorList) + 1):
			tempNotes = self.currentNotes[:]
			self.double_transpose(*generatingVL)
			vlPaths.append([self.currentNotes[j] - tempNotes[j] for j in range(len(self.currentNotes))])
		
		"""apply the circle-of-fifths voice motions to the original chord"""
		self.currentNotes = tempPCs
		self.lastNotes = tempPCs
		
		self.latticeApplications = len(vectorList)
		totalPaths = [sum([vlPaths[j][i] for j in vectorList]) for i in range(len(self.currentNotes))]
		self.lastVoiceleading = {(self.currentNotes[i] % abs(self.modulus)):totalPaths[i] for i in range(len(self.currentNotes))}
		self.suppressChangeFunc = changeFuncStatus
		self.apply_VL_dict(self.lastVoiceleading)
	
	def undo_lattice(self, suppressChange = True):
		
		"""
		undoes the last lattice voice leading, then applies the basic VL to move to the chord with the same sum
			aka "move to the central spine of the lattice
		"""
		if not self.latticeApplications:
			return
		oldChangeFunc = self.suppressChangeFunc
		if suppressChange:
			self.suppressChangeFunc = True
		self.undo()										# return to the place we were before the lattice moves
		for i in range(self.latticeApplications):
			self.double_transpose(*self.latticeVL)		# generating VLs equal to the number of lattice moves
		self.latticeApplications = 0					# reset to indicate we are at the central spine
		if suppressChange:
			self.suppressChangeFunc = oldChangeFunc
	
	def shift_lattice(self, bvlMoves = 0):
		"""moves along the latice by bvlMoves, preserving lattice position"""
		self.undo_lattice()
		self.double_transpose(self.basicVL[0]*bvlMoves, self.basicVL[1]*bvlMoves)
		self.redo_lattice()
	
	def redo_lattice(self):
		if self.lastLatticeParameters:
			self.lattice(*self.lastLatticeParameters)
			
	def find_transformation(self, startPitches = None, endPitches = None):
		
		if not startPitches and not endPitches and self.lastPermutation:
			"""TODO: need to add code here!!!!"""
			pass
			
		if not startPitches:
			startPitches = self.lastNotes
		
		if not endPitches:
			endPitches = self.currentNotes
	
	"""========= find the permutation and transposition resulting from any voice leading ========="""			
	def find_transposition_and_permutation(self, startPitches = None, endPitches = None, invert = False):
		"""
		add code to deal with neighbors, inversions, and pedals
		
		check inversion sign
		
		"""
		
		if not startPitches:
			startPitches = self.lastNotes
		
		if not endPitches:
			endPitches = self.currentNotes
		
		if invert:
			sign = -1
		else:
			sign = 1
			
		l = len(startPitches)
		theInts = [[(endPitches[i] - sign*startPitches[j]) % self.modulus for j in range(l)] for i in range(l)] 			# list of ints from note i of chord 1 to all notes of chord 2
		targetInts = list(set(theInts[0]).intersection(*theInts[1:]))														# what ints are common to all?										
		
		if not targetInts: 
			return None, None
			
		targetInt = targetInts[0]
		thePerm = []
		theTransps = []
		
		for n in startPitches:
			for i, n2 in enumerate(endPitches):
				theInt = n2 - sign*n
				if theInt % self.modulus == targetInt:
					thePerm.append(i)
					theTransps.append(theInt)
					break
					
		if len(list(set(thePerm))) != l:
			return None, None
			
		permDict = {i:thePerm[i] for i in range(len(thePerm))}
		transpDict = {i:theTransps[i] for i in range(len(theTransps))}
		
		return permDict, transpDict
	
	"""========= SUPPORT NONHARMONIC TONES ========="""
	"""TODO: get code from Arpeggiation"""
			
	def __iter__(self):
		self.iterIndex = 0
		return self
	
	def __next__(self):
		if self.iterIndex < len(self.currentNotes):
			self.iterIndex += 1
			return self.currentNotes[self.iterIndex - 1]
		raise StopIteration
	
	def __len__(self):
		return len(self.currentNotes)
		
	def rectify_multiset(self, pitchList):
		
		fullStructure = copy.deepcopy(pitchList)
		pcCounts = {}
		
		foundMultiset = False
		
		for i, n in enumerate(pitchList):
			
			if type(n) is not list:
				PC = n % self.modulus
				if type(PC) is float:
					PC = round(PC, 2)
					n = round(n, 2)
					
				count = pcCounts.get(PC, -1) + 1
				pcCounts[PC] = count
				
				if count:
					fullStructure[i] = round(n + .0001*count, 4)				# turn duplicate into fractionally displaced PC
					foundMultiset = True
		
		if foundMultiset:
			return fullStructure
		
		return False

	"""======== TEXT TRANSFORMATIONS ====================="""
	
	def transform(self, tString, getCallChain = False):
		
		"""
		
		Possible commands:
		
			T - transpose along scale
			t - transpose along chord
				TODO: T* and t* for dualist transposition (currently S) [or choose some other symbol for dualist transformations]
		
			C - crossing or voice exchange
				C(x, y) - exchange notes x/y (order matters: the first goes up to the second, the second goes down to the first)
				C3		- crossing number 3
				TODO: C* dualist crossings
		
			I - traditional inversion within scale
			i - traditional inversion within chord
		
			J - NR inversion (combines I and i)
			 	J(A, B) - NR inversion around AB
				Jx		- NR inversion number x (nondualist)
				TODO: J* - NR inversion number x (dualist)??
		
			Kx - NR inversion number x (nondualist) TODO: rename
		
			S - dualist transposition TODO: rename
		
			L - motion along the 'voice leading lattice' (Tymoczko "Scale Networks")
			
			P(3, -1, 2)[1, 2, 0] - voice 0 is transposed by 3 and moves to voice 1, voice 1 is transposed by -1 and moves to voice 2, etc.
			Q(3, -1, 2)[1, 2, 0] - voice 0 moves by melodic interval 3 and then its material moves to voice 1, etc.
			
			M(A: -1, B: -2, C: 3) - simple voice leading; can also use VL(A: -1, B: -2, C: 3)
			V(A: -1, B: -2, C: 3) - transpositional VL schema
			W(A: -1, B: -2, C: 3) - transpositional contrapuntal interchange
			X(A: -1, B: -2, C: 3) - inversional contrapuntal interchange
			TODO: 
				need inversional contrapuntal VL schema 
				all contrapuntal interchanges should have a counter (maybe they do?)
				once defined you should be able to reuse contrapuntal interchanges without relabeling them
					in other words, having done c.transform('X1(0:7, 4:-1, 7:-7)') you should be able to do c.transform('X1')
		
		TODO: 
			Need dualist crossings
			Name the dualist transforms systematically?  J* and C* for dualist?
		
		"""
		callChain = []
		
		"""Text formatting:
		
			1. remove duplicate spaces
			2. convert to uppercase except for those transforms that are case sensitive
			3. eliminate spaces after colons (for the VL schemas)
		
		"""
		
		# NB: there is some weirdness with the cases here, because i and t end up lowercase. check out CoreVoiceleading.transformSubstitutions
		tString = ' '.join(tString.split(' '))
		tString = tString.replace('t', 'TTTT').replace('i', 'IIII').upper().replace('TTTT', 't').replace('IIII', 'i').replace(': ', ':').replace(' ', ',')
		for k, v in CoreVoiceleading.transformSubstitutions.items():
			tString = tString.replace(k, v)
		
		tString = ','.join([x for x in tString.split(',') if x])							# split works differently for spaces and commas!
		
		commandList = parse_parentheses(list(tString), joinSymbol = '')						# group by parentheses
		
		commandLetters = 'CIiJKLMPQSTtVWX'													# TODO: check these
				
		curStr = ''
		
		for c in commandList + ['z']:														# add a string terminator
			
			if c[0] == 'z' or c[0] in commandLetters:										# we have reached a new command, process the last
				
				out = self.parse_individual_text_transformation(curStr)
				
				if out:
					callChain.append(out)
					
				curStr = c
			else:
				curStr += c
		
		if getCallChain:
			return callChain
		
		return self.apply_call_chain(callChain)
	
	def apply_call_chain(self, callChain):
		
		if not callChain:
			return False
		
		lastNotes = self.currentNotes[:]						# store the initial state
		self.lastTransposition = 0
		self.lastPermutation = []
		self.lastPaths = []
		self.lastLabel = None
			
		for func, args, kwargs in callChain:
			func(*args, **kwargs)
			
		self.outPCs = self.currentNotes
		self.lastNotes = lastNotes
		
		return self.currentNotes
	
	def parse_individual_text_transformation(self, curStr):
		
		"""
		
		TODO:
			change attribute names, and make each transformation return the main object, so you can type
				myObject.T(4).t(-1)
		
		"""
		
		if len(curStr) <= 1: return False
		
		tType = curStr[0]
		tLabel = curStr[1:]
	
		
		if tType == 'T':
			tLabel = ''.join([x for x in tLabel if (x.isdigit() or x == '-')])
			if tLabel:
				data = int(tLabel)
				self.lastTransposition += data
				return [self.double_transpose, [data, 0], {}]
				
		elif tType == 't':
			tLabel = ''.join([x for x in tLabel if (x.isdigit() or x == '-')])
			if tLabel:
				return [self.double_transpose, [0, int(tLabel)], {}]
		
		elif tType == 'S':
			extraRotations = tLabel.count('.')
			tLabel = ''.join([x for x in tLabel if (x.isdigit() or x == '-')])
			if tLabel:
				return [self.dualist_transposition, [int(tLabel), extraRotations], {}]
		
		elif tType == 'V':
			pairs = self.get_schema_from_text(tLabel)
			if pairs:
				return [self.transpositional_VL_schema, [pairs], {}]
		
		elif tType == 'W':
			pairs = self.get_schema_from_text(tLabel)
			if pairs:
				listOfSchemas = [pairs] + [self.reverse_schema(pairs)]
				#print('here', listOfSchemas)
				return [self.multischema, [listOfSchemas], {}]
		
		elif tType == 'X':
			pairs = self.get_schema_from_text(tLabel)
			if pairs:
				iPairs = self.invert_schema(pairs)
				listOfSchemas = [pairs] + [self.reverse_schema(pairs)] + [iPairs] + [self.reverse_schema(iPairs)]
				#print('here', listOfSchemas)
				return [self.multischema, [listOfSchemas], {}]
				
		elif tType == 'M':
			pairs = self.get_schema_from_text(tLabel)
			if pairs:
				return [self.voicelead, [pairs], {}]
		
		elif tType == 'C':															# uses a different internal label from the text command, which is C.
			extraRotations = tLabel.count('.')
			if tLabel.count('('):
				PCs = self.find_inversion_from_text(tLabel)							
				if self.container:
					PCs = [self.container.scalarPC(x) for x in PCs]
				PCs = [x % self.modulus for x in PCs]
				if PCs[1] <= PCs[0]: PCs[1] += self.modulus
				if len(PCs) == 2:
					pairs = [[PCs[0], PCs[1] - PCs[0]], [PCs[1], PCs[0] - PCs[1]]]
					return [self.voicelead, [pairs], {}]
			else:
				tLabel = ''.join([x for x in tLabel if (x.isdigit() or x == '-')])
				if tLabel:
					return [self.apply_crossing, [int(tLabel), extraRotations], {}]
					
		elif tType == 'I' or tType == 'i':
			if tType == 'i':
				inversionFunction = self.invert_within_chord
				indexNumberContainer = self.scale
			else:
				inversionFunction = self.invert
				indexNumberContainer = self.container
			if tLabel.count('('):
				pitches = [parse_lettername(x) for x in tLabel.split(',')]
				if len(pitches) == 1:
					pitches.append(pitches[0])						
				if self.container:
					pitches = [self.container.scale_degree(x) for x in pitches]
				if len(pitches) >= 2:
					return [inversionFunction, [pitches[0], pitches[1]], {}]
			else:
				"""
				index numbers are always specified chromatically, in keeping with the principle of avoiding scale degree numbers
				the way they work is: they are divided by two to get a fixed point, and then quantized to the nearest upper and lower collection tones (which can be identical)
				this means that they might work differently in different collections: i122 is inversion around C4 E4 within the C major triad but C D within the C major scale
				"""
				chromaticIndexNumber = parse_number(tLabel)
				if indexNumberContainer:
					pitches = self.get_pair_from_index_number(chromaticIndexNumber, indexNumberContainer)
				else:
					pitches = [chromaticIndexNumber, 0]
				return [inversionFunction, [pitches[0], pitches[1]], {}]
				
		elif tType == 'J':
			if tLabel.count('('):
				PCs = self.find_inversion_from_text(tLabel)							
				if self.container:
					PCs = [self.container.scalarPC(x) for x in thePitches]
				PCs = [x % self.modulus for x in PCs]
				if len(PCs) == 2:
					return [self.NR_inversion, [PCs], {}]
			else:
				extraRotations = tLabel.count('.')
				tLabel = ''.join([x for x in tLabel if (x.isdigit() or x == '-')])
				if tLabel:
					return [self.NR_inversion_by_number, [int(tLabel), extraRotations, True], {}]
		
		elif tType == 'K':
			extraRotations = tLabel.count('.')
			tLabel = ''.join([x for x in tLabel if (x.isdigit() or x == '-')])
			if tLabel:
				return [self.NR_inversion_by_number, [int(tLabel), extraRotations, False], {}]
		
		elif tType == 'L':
			if tLabel.count('*'):
				return [self.undo_lattice, [], {}]
			latticeMoves = [parse_number(x) for x in tLabel.split()]
			return [self.lattice, latticeMoves, {}]
		
		elif tType == 'P':
			self.lastLabel = curStr
			repeatVL = self.permDict.get(self.lastLabel, False)
			if repeatVL:
				paths, permutation = self.get_perturbations(repeatVL)
			else:
				paths, permutation = self.get_perturbations(tLabel)
			if paths:
				return [self.apply_perturbation, [paths, permutation], {}]
			
		elif tType == 'Q':
			self.lastLabel = curStr
			repeatVL = self.permDict.get(self.lastLabel, False)
			if repeatVL:
				paths, permutation = self.get_perturbations(repeatVL)
			else:
				paths, permutation = self.get_perturbations(tLabel)
			if paths:
				return [self.apply_perturbation_melodic, [paths, permutation], {}]
	
	"""TODO: complete this list, so that transformations can be accessed wtihout text translation, in a consistent way"""
	def T(self, param = 0):
		self.double_transpose(param, 0)
		return self								# needed to chain transformations together, as in a.T(4).t(-1)
	
	def t(self, param = 0):
		self.double_transpose(0, param)
		return self
		
	def S(self, *args):
		self.dualist_transposition(*args)
		return self
	
	def V(self, *args):
		self.transpositional_VL_schema(*args)
		return self
	
	def C(self, *args):
		self.apply_crossing(*args)
		return self
		
	def apply_perturbation(self, paths, permutation): 
		self.lastPaths = paths
		
		outPCs = [None] * len(self.currentNotes)
		
		if type(permutation) is str and permutation.count('CE'):
			nextPaths = [-x for x in paths]
			self.permDict[self.lastLabel] = 'P(' + str(nextPaths)[1:-1] + ')' + str(permutation)
			permutation = list(range(len(paths)))
		
		tDict = {(self.currentNotes[i] % self.modulus):paths[i] for i in range(len(paths))}
		
		self.lastPermutation = permutation	
		
		"""TODO: figure out if I meant anything when I wrote this? it doesn't make sense"""
		"""we might want to save the lastTransposition at some point, but that should happen elsewhere, not here"""
		"""we might want to save these voice leadings but we should give them a different name"""
		#self.lastTransposition = tDict
		
		self.apply_VL_dict(tDict)
	
	def apply_perturbation_melodic(self, melodyPaths, permutation):
		
		paths = [None] * len(self.currentNotes)
		
		for i, melodicInterval in enumerate(melodyPaths):
			sourceVoiceNumber = permutation.index(i)
			newPath = self.currentNotes[i] + melodicInterval - self.currentNotes[sourceVoiceNumber]
			paths[sourceVoiceNumber] = newPath
		
		self.permDict[self.lastLabel] = 'P(' + str(paths)[1:-1] + ')' + str(permutation)			# these don't repeat in an intuitive way; they need to be converted to regular permutations
		self.apply_perturbation(paths, permutation)
		
	def NR_inversion(self, iList):
		self.apply_VL_dict(NR_inversion_dict(self.sortedPCs, iList, self.modulus))
	
	def NR_inversion_by_number(self, inversionNumber = 0, extraRotations = 0, dualistLabels = False):
		self.apply_VL_dict(NR_inversion_by_number_dict(self.sortedPCs, inversionNumber, extraRotations, dualistLabels, self.modulus))
	
	def get_pair_from_index_number(self, chromaticIndexNumber, collection):
		fixedPoint = chromaticIndexNumber/2
		upperNote = collection.scale_degree(fixedPoint, direction = 1, recurse = True)
		lowerNote = collection.scale_degree(fixedPoint, direction = -1, recurse = True)
		if collection.container:
			lowerNote = collection.__getitem__(lowerNote, recurse = False)
			upperNote = collection.__getitem__(upperNote, recurse = False)
		return lowerNote, upperNote
	
	def invert(self, note1, note2):
		"""traditional inversion"""
		indexNumber = note1 + note2
		invertedNotes = [indexNumber - x for x in self.currentNotes]
		tDict = {self.currentNotes[i]%self.modulus:invertedNotes[i]-self.currentNotes[i] for i in range(len(self.currentNotes))}
		self.apply_VL_dict(tDict)
	
	def invert_within_chord(self, note1, note2):
		originalSpacing = self.get_spacing([note1, note2] + self.currentNotes, note1)
		spacing = [-x for x in originalSpacing]
		invertedNotes = self.apply_spacing(spacing, note2)[2:]
		tDict = {self.currentNotes[i]%self.modulus:invertedNotes[i]-self.currentNotes[i] for i in range(len(self.currentNotes))}
		self.apply_VL_dict(tDict)
	
	def make_Tt_dict(self, bigT = 0, littleT = 0):
		tDict = {}
		pSort = self.sortedPCs + [x + self.modulus for x in self.sortedPCs]
		octaves = int((littleT + EPSILON) / len(self.currentNotes))
		if littleT < 0:
			octaves = octaves - 1
		bigT += self.modulus * octaves
		steps = littleT % len(self.currentNotes)
		for i, pc in enumerate(self.sortedPCs):
			path = pSort[i+steps] - pc + bigT
			tDict[pc] = path
		return tDict
	
	def double_transpose(self, bigT = 0, littleT = 0):
		self.apply_VL_dict(self.make_Tt_dict(bigT, littleT))
	
	def apply_crossing(self, crossingNumber = 0, extraRotations = 0):
		tDict = {}
		gnf = geometrical_normal_form_local(self.sortedPCs, invert = False, extraRotations = extraRotations, modulus = self.modulus)
		transp = None
		for i in range(self.modulus):
			if all([(x + i) % self.modulus in self.sortedPCs for x in gnf]):
				transp = i
				break
		if transp is None: return {}
		tDict = {}
		v1 = (-crossingNumber) % len(self.sortedPCs)
		v2 = (v1 + 1) % len(self.sortedPCs)
		note1 = (gnf[v1] + transp) % self.modulus
		note2 = (gnf[v2] + transp) % self.modulus
		if note2 < note1:
			note2 += self.modulus
		path1 = note2 - note1
		path2 = -path1
		tDict = {note1: path1, (note2 % self.modulus): path2}
		self.apply_VL_dict(tDict)
		
	def dualist_transposition(self, transposition = 0, forceInversions = False):
		"""inversion by number, either dualist or nondualist"""
		
		gnf = geometrical_normal_form_local(self.sortedPCs, invert = True, modulus = self.modulus)
		cardinality = len(gnf)
		
		transp = None
		
		for i in range(self.modulus):
			if all([(x + i) % self.modulus in self.sortedPCs for x in gnf]):
				transp = i
				break
				
		if (not forceInversions) and (transp is not None):
			self.apply_VL_dict({}, globalTransposition = transposition)
			return
		
		index = None
		for i in range(self.modulus):
			if all([(i - x) % self.modulus in self.sortedPCs for x in gnf]):
				index = i
				break
				
		if index is None:
			return False
			
		self.apply_VL_dict({}, globalTransposition = -transposition)
	
	def get_perturbations(self, inputSym):
		"""raw permutations ordered by voices"""
		sUp = inputSym.upper().replace(';', ',').replace('(', '').replace(')', '').replace(' ', '')
		while sUp.count(',,'):
			sUp = sUp.replace(',,', ',')
		permutation = list(range(len(self.currentNotes)))
		if sUp.count('['):
			sUpSplit = sUp.split('[')
			sUp = sUpSplit[0]
			permString = sUpSplit[1]
			if permString.count('CE') or permString.count('<->'):					# contrapuntal interchange, back and forth
				permutation = '[CE]'
			elif permString.count(','):												# specific permutation to apply
				permutation = [parse_number(x) for x in permString.split(',')]
		sSplit = sUp.split(',')
		paths = [parse_number(x) for x in sSplit if x]
		return paths, permutation
		
	def get_schema_from_text(self, inputSym, useMidi = True, fill = False):
		sUp = inputSym.upper().replace(';', ',').replace('(', '').replace(')', '').replace(' ', '').replace('V', '')
		out = {}
		for s in sUp.split(','):
			if not s.count(':'): continue
			s = s.split(':')
			pc = parse_lettername(s[0]) % 12
			path = parse_number(''.join([x for x in s[1] if x in '-.' or x.isdigit()]))
			if useMidi and self.container:
				pc = self.container.scalarPC(pc)		
			out[pc] = path
		"""for pc in self.PCs:
			if pc not in out:
				if useMidi and self.container:
					pc = self.container.scalarPC(pc)
				out[pc] = 0"""
		if fill:
			for p in self.PCs:
				if p not in out: out[p] = 0
		return list(out.items())
	
	def find_inversion_from_text(self, tLabel):
		pcPairs = tLabel.replace('(', '').replace(')', '').split(',')
		if len(pcPairs) == 1:
			pcPairs.append(pcPairs[0])
		if len(pcPairs) > 2:
			pcPairs = pcPairs[:2]
		if len(pcPairs) != 2: 
			return False
		pcPairs = [parse_lettername(x) for x in pcPairs]
		return pcPairs
	
	@staticmethod
	def get_transpositions(n, o):
	
		"""
		returns the basic voice leading transpositions for an n-note chord in an o-note scale expressed as [T = scalesteps, t = chordsteps]
			when n divides o, we have a sum-0 voice leading with t = 1
			when n and o are relatively prime it returns the sum -1 voice leading
			when n and o share a common factor it returns a pair of voice leadings, the sum-0 voice leading and the smallest descending voice leading
		
		"""
		
		if o == 1 or n == 1:
			return [[1, 0]]
			
		t = int(o/n)
		if t*n == o:
			return [[-t, 1]]
		for i in range(1, o):
			t = n * i
			if t % o == 1:
				return [[o-i, -((-int(t/o)) % n)]]
			elif t % o == 0:
				break
		multiplicity = int(o/i)
		vl1 = [i, -int(t/o)]
		for i in range(o-1, 1, -1):
			t = n * i
			if t % o == multiplicity:
				return [vl1, [o-i, -((-int(t/o)) % n)]]
		print("ARCA ERROR: can't find basicVL: ", o, n)
		return [[1, 0]]

"""silly little class to allow individual indexing of the repeat function call chain while still using a setter to perform translation"""		
class RepeatFunctionWrapper():
	
	def __init__(self, container, initialValue = None):
		self.container = container
		if initialValue is not None:
			self.set_call_chain(initialValue)
	
	def set_call_chain(self, theList):
		self.container.repeatCallChain = theList
	
	def __getitem__(self, k):
		return self.container.repeatCallChain[k]
	
	def __setitem__(self, k, val):
		self.container.repeatCallChain[k] = self.container.translate_repeat_function(val)
	
	def __str__(self):
		return str(self.container.repeatCallChain)
	
	def __iter__(self):
		self.iterIndex = 0
		return self
	
	def __next__(self):
		if self.iterIndex < len(self.container.repeatCallChain):
			self.iterIndex += 1
			return self.container.repeatCallChain[self.iterIndex - 1]
		raise StopIteration
	
	def __len__(self):
		return len(self.container.repeatCallChain)
		
def progression_to_round(*args, **kwargs):
	
	print("Progression_to_round has been moved to Scale.  Use\n\nbasicRound = arca.Scale(chordList = ['C4 Eb4 F4 Bb4', 'Bb3 D4 G4 A4', 'Bb3 C4 F4 Ab4', 'Bb3 C4 Eb4 F4'])\n\ninstead of"+
	"\n\nbasicRound = arca.progression_to_round([[0, 3, 5, 10], [-2, 2, 7, 9], [-2, 0, 5, 8], [-2, 0, 3, 5]])")

class Scale(CoreVoiceleading):
	
	defaultScales = {7: [60, 62, 64, 65, 67, 69, 71]}
	
	defaults = {		'container': None, 
						'modulus': None, 
						'keyOffset': 60, 
						'startDegree': 60, 
						'progression': [], 
						'useMidi': True, 
						'sort': True,					# sort the scale's notes (for useMidi)
						'changeFunc': None,
					}
	
	"""
	
	A Nestable Scale
	
		NB: scale need not be sorted, though fractional scale degrees will need a sorted scale
	
	The core routines are
	
		__getitem__ :	gives the scale degree of the containing scale corresponding to an input scale degree; 
						with recurse = True, chases through the hierarchy to get a MIDI note
	
		scale_degree :	the reverse function; accepts a scale degree or MIDI note (recurse=True) and delivers the closest 
						input scale degree

	When you initialize, you choose keyOffset and scartDegree.  I have found it is considerably more intuitive to use numbers keyed to the middle C octave, 
		rather than PCs, so keyOffset = 60 and scaleDegree = 60.
	
		For input pitches [a, b, c, ...] this will map scale degrees [startDegree, startDegree + 1, startDegree + 2, ...] to [a + keyOffset, b + keyOFfset, c + keyOffset ...]
	
	modulus is normally the modulus of the containing scale but sometimes you might want to change it:
		modulus = -12 will create a scale where ascending steps map to descending intervals (inverting)
		modulus = 11 can create a nonoctave repeating scale in the 12-tone chromatic scale
	
	I believe this works with floating point numbers, but that has not been carefully tested
					
	Can initialize with either chordList = [...] for a round/repeating progression, or sequence = [...] for a sequence.  These are not default args.
	
	"""
	
	def __init__(self, *args, **kwargs):
		self.currentNotes = []
		self.initialize_core()
		self.change_scale(*args, **kwargs)
		self.get_basic_VL()
		self.initialNotes = self.currentNotes[:]
		self.scale = self											# for consistency with the Chord class 
		#self.__getitem__ = self.__getitemBASE__
	
	def copy(self, scaleObj):
		self.change_scale(scaleObj.PCs, scaleObj.container, scaleObj.modulus, 0, scaleObj.startDegree, scaleObj.progression)
		
	def update_scale(self, *args, **kwargs):
		"""meant as a lighter and more flexible version of the routine below"""
		
		for attribute, value in kwargs.items():
			if hasattr(self, attribute):
				setattr(self, attribute, value)
				
		if hasattr(self, 'chordList') or hasattr(self, 'sequence'):
			if hasattr(self, 'sequence'):
				startObject, prog = self.chords_to_transformations(self.sequence, schema = True)
			else:
				startObject, prog = self.chords_to_transformations(self.chordList, schema = False)
			startObject = sorted([x % self.modulus for x in startObject])
			args = [startObject]
			self.progression = prog
		
		if args:
			theScale = self.check_for_multiset(args[0])
			self.len = len(theScale)
			self.lastVoiceleading = {}
			self.latticeApplications = 0														# for undoing VL-lattice moves
			self.lastLatticeParameters = []														# for redoing lattice moves
			self.currentNotes = [x + self.keyOffset for x in theScale]
			self.SDs = [self.startDegree + i for i in range(self.len)] 
			self.modSDs = [x % self.len for x in self.SDs]
			self.update()
			self.get_basic_VL()
	
	def check_for_multiset(self, theScale):
		temp = self.rectify_multiset(theScale)
		if temp:
			self.output_note = self.output_note_multiset
			theScale = temp
		else:
			self.output_note = self.output_note_regular
		return theScale
		
	def spell(self):
		return ' '.join(spell_MIDI(self.currentNotes, octaves = False))
		
	def change_scale(self, *args, **kwargs):
		
		self.set_default_arguments(kwargs)
		
		if self.modulus is None:
			self.modulus = len(self.container) if self.container else 12
			
		if hasattr(self, 'chordList') or hasattr(self, 'sequence'):
			if hasattr(self, 'sequence'):
				startObject, prog = self.chords_to_transformations(self.sequence, schema = True)
			else:
				startObject, prog = self.chords_to_transformations(self.chordList, schema = False)
			startObject = sorted([x % self.modulus for x in startObject])
			args = [startObject]
			self.progression = prog
		
		if args:
			theScale = args[0]
			
		self.hierarchy = self.get_hierarchy()
		
		self.permDict = {}
		
		if type(theScale) is int:
			theScale = Scale.defaultScales.get(theScale, [int(i*self.modulus/theScale) for i in range(theScale)])	# maximally even n-note scale
		elif type(theScale) is str:
			theScale = [parse_lettername(x) for x in theScale.split()]
			
		if self.useMidi:	
			if self.container:
				theScale = [self.container.scalarPC(x) for x in theScale]
			else:
				theScale = [x % 12 for x in theScale]	
				
		if self.sort:
			theScale = sorted(theScale)
		
		theScale = self.check_for_multiset(theScale)
		
		self.len = len(theScale)
		self.lastVoiceleading = {}
		self.latticeApplications = 0														# for undoing VL-lattice moves
		self.lastLatticeParameters = []														# for redoing lattice moves
		
		self.currentNotes = [x + self.keyOffset for x in theScale]
		self.SDs = [x + self.startDegree for x in range(self.len)] 
		self.modSDs = [x % self.len for x in self.SDs]
		
		self.update()
		self.get_basic_VL()
		
	def change_scale_with_MIDI(self, theNotes, **kwargs):
		PCs = sorted(set([x % 12 for x in theNotes]))
		if self.container:
			PCs = [self.PC(x) for x in PCs]
		self.change_scale(PCs, container = self.container, **kwargs)
	
	"""
	Converting between MIDI and scale degrees, and pitch and pitch class, is a persistent source of confusion.
		
		MIDI -> chromatic PC				PC -> x % 12
		MIDI -> scalar PC 					((x % 12) + obj.keyOffset) % obj.modulus
		MIDI -> scale degree 				obj.scale_degree(x, recurse = True)		# recurse is true by default
				``
		scale degree -> MIDI				obj.__getitem__(x, recurse = True)		# recurse is true by default
		scale degree -> chromatic PC		obj.__getitem__(x, recurse = True) % 12	# recurse is true by default
		scale degree -> scalar PC			x % obj.modulus
		
	"""
	
	def scalarPC(self, midiNote):
		return self.scale_degree((midiNote % 12) + self.keyOffset) % self.modulus			# this is very confusing, but needed
	
	def scalarPC_internal(self, scaleDegree):
		return (scaleDegree - self.keyOffset) % self.len									# version of this for internal scale degrees
	
	def quantize(self, n, direction = 0, tieBreak = 1, extraTransposition = 0):
		return self.__getitem__(self.scale_degree(n, direction, tieBreak, recurse = True) + extraTransposition)
	
	def update(self):
		self.PCs = [x % abs(self.modulus) for x in self.currentNotes]						# voicelead can produce negative numbers or numbers > self.modulus
		self.sortedPCs = sorted(self.PCs)													# needed for quantizing
		self.wrappedSortedPCs = self.sortedPCs + [self.sortedPCs[0] + abs(self.modulus)]	# useful for scale-degree calculations
		self.doublePCs = self.sortedPCs + [x + abs(self.modulus) for x in self.sortedPCs]	# needed for big T little t
		self.sortedSDs = [self.SDs[self.PCs.index(x)] for x in self.sortedPCs]
	
	"""simple routine to reset the notes of the secale without changing much else"""
	def set_notes(self, theScale):
		self.currentNotes = [x + self.keyOffset for x in theScale]
		self.update()
		
	def double_transpose(self, T = 0, t = 0):
		"""make a VLdict for transposition along intrinsic and extrinsic scale"""
		if t < 0:
			oShift = int((t + EPSILON)/self.len) - 1
			T = T + oShift*self.modulus
			t = t % self.len
		if t > self.len:
			oShift = int(t/self.len)
			T = T + oShift*self.modulus
			t = t % self.len
		if t == 0:
			self.apply_VL_dict({(self.sortedPCs[i]%self.modulus):T for i in range(self.len)})
		else:
			self.apply_VL_dict({(self.sortedPCs[i]%self.modulus):(T + self.doublePCs[i+t] - self.sortedPCs[i]) for i in range(self.len)})
	
	def transpose_note(self, n = 60, transp = 0):
		sd = self.scale_degree(n)
		return self.__getitem__(sd+transp)
	
	def make_voicing(self, intervalList = [], **kwargs):
		
		stepList = []
		theSum = 0
		
		for i in intervalList:
			theSum += i
			stepList.append(theSum)
		
		return self.space_chord(stepList = stepList, **kwargs)
	
	def space_chord(self, stepSize = 1, adjustment = 1, separatingSteps = 1, topNote = 0, bottomNote = 0, stepList = [], level = None):
		
		"""
		
		spaces a chord according to stepSize (1 = close position, 2 = open position, etc.)
		
		when the stepSize divides the chord evenly, you have multiple separate cycles; adjustment adjusts the interval between the cycles
		when separatingSteps is 0 adjustment applies directly to the top note, e.g. C4-G4-F#4-B4 for adjustment -1 and spacing 2.
		when separatingSteps is 1, adjustment applies to the top note + one spacing step, e.g. C4-G4-B4-F#5 for adjustment -1 and spacing 2.
		
		"""
				
		if stepList:
			stepList = [0] + stepList					# I feel weird about this!
			
		else:
			if stepSize > self.len/2:
				originalCycle = self.len/(self.len-stepSize)
			else:
				originalCycle = self.len/stepSize
			cycleSize = int(originalCycle)
			if cycleSize == originalCycle:
				numCycles = int(self.len/cycleSize)
				stepList = [i * stepSize for i in range(cycleSize)]
				for i in range(numCycles - 1):
					startNote = stepList[-1] + adjustment + (stepSize*separatingSteps)
					stepList += [(i * stepSize)+startNote for i in range(cycleSize)]
			else:
				stepList = [i*stepSize for i in range(self.len)]
		if topNote:
			adjustment = topNote - stepList[-1]
		elif bottomNote:
			adjustment = bottomNote - stepList[0]
		else:
			return [self.__getitem__(x, depth = level) for x in stepList]
		return [self.__getitem__(x + adjustment, depth = level) for x in stepList]
	
	def get_repeat_function_permutation(self):
		
		"""
		how does the total collection of repeat functions transform the chord's voices?
		does not yet work with custom referenceList
		
		"""
		
		perm = list(range(self.len))
		
		for fList in self.progression:
			if type(fList) is list:
				if type(fList[0]) is str:
					fList = [fList]
				for f in fList:
					if f[0] == 'Tt':
						sTransp = f[1][1]
						perm = [(x + sTransp) % self.len for x in perm]
					elif f[0] == 'crossing':
						voices = f[1]
						v1 = voices[0]
						if len(voices) == 2:
							v2 = voices[1]
						else:
							v2 = (v1 + 1) % self.len
						perm = [v1 if x == v2 else (v2 if x == v1 else x) for x in perm]
		return perm
	
	def get_hierarchy(self):
		"""get the hierarchy of scales below the current one"""
		containers = [self]
		while containers[-1].container:
			containers.append(containers[-1].container)
		return containers
	
	def get_note_in_hierarchy(self, note, level = None, depth = None):
		if level is None or level >= len(self.hierarchy):
			return note
		collection = self.hierarchy[level % len(self.hierarchy)]
		return collection.__getitem__(note, depth = depth)
	
	def get_interval_in_hierarchy(self, note, targetNote, level = None, checkErrors = True):			# specify notes with MIDI; this is for nonharmonic tonesf
		if level is None or level >= len(self.hierarchy):
			scalarInterval = (targetNote - note) % 12
			l = 12
		else:
			collection = self.hierarchy[level % len(self.hierarchy)]
			l = len(collection)
			if checkErrors:
				n = collection.notes
				if note % 12 not in n: 
					print("Pitch class", note %12, "not in scale", n)
					return
				elif targetNote % 12 not in n: 
					print("Pitch class", targetNote %12, "not in scale", n)
					return
			scalarInterval = (collection.scale_degree(targetNote) - collection.scale_degree(note)) % l
		if scalarInterval > l/2: return scalarInterval - l
		return scalarInterval
	
	def has_PCs(self, pcList):
		rc = self.raw_content()
		return all([x % 12 in rc for x in pcList])
	
	def make_quantization_VL(self, shiftDir = 0, tieBreak = 0):
		vlDict = {}
		for i in range(0, 12):
			vlDict[i] = self.__getitem__(self.scale_degree(i + 60, shiftDir = shiftDir, tieBreak = tieBreak), recurse = True) - 60 - i
		return vlDict
	
	def __contains__(self, item):
		return item % abs(self.modulus) in self.sortedPCs
	
	@property
	def notes(self):
		return self.raw_content()
	
	def __getitem__(self, key, recurse = True, depth = None):
		
		"""
		
		with recurse or depth is None: get the midivalue corresponding to a scaledegree
		without recurse: get the scale degree of the containing scale corresponding to a scaledegree
		with depth: go n levels deep in the hierarchy
			0 -- surface level of hierarchy (self, the current scale)
			1, 2, ... for the first, second, etc, scale in the hierarchy
			-1, -2, -3 for the final, second to last, etc., scale in the hierarchy
		
		"""
		
		if isinstance(key, slice):
			start = key.start
			end = key.stop if (key.stop is not None) else start + len(self.currentNotes)
			step = key.step if (key.step is not None) else 1
			return [self.__getitem__(x, recurse = recurse, depth = depth) for x in range(start, end, step)]
		
		if depth == 0: return self.output_note(key)
		
		if recurse:
			depth = None
			
		if depth:
			depth = depth % len(self.hierarchy)
		
		if type(key) is int:
			index = self.modSDs.index(key % self.len)
			PC = self.currentNotes[index]
			SD = self.SDs[index]
			neededOctaves = int((key - SD)/self.len)
			output = PC + neededOctaves * self.modulus
			if depth:
				if depth > 1 and self.container:
					return self.output_note(self.container.__getitem__(output, recurse = False, depth = depth - 1))
				else:
					return self.output_note(output)
			elif recurse and self.container:
				return self.output_note(self.container.__getitem__(output, recurse = True))
			else:
				return self.output_note(output)
		
		k = int(key)
		frac = key - k
		kVal1 = self.__getitem__(k, recurse = recurse, depth = depth)
		kVal2 = self.__getitem__(k + 1, recurse = recurse, depth = depth)
		return self.output_note(kVal1 + (kVal2 - kVal1)*frac)
		
	def scale_degree(self, target, direction = 1, tieBreak = 0, recurse = True):
		"""
		
		with recurse: find the scale degree closest to a midi note 
		without recurse: find the scale degree closest to a containing scale degree
		
			direction is either > 0 for up, < 0 for down, 0 for closest, or None for don't quantize (fractional value)
		
		Since a scale can be arbitrarily ordered, the non-quantizing option is somewhat tricky:
			(1) if the target lies between two adjacent scale tones a fractional value is returned
			(2) otherwise it returns a list of the neearest scale tones [lowerTone, upperTone]
		
		tieBreak controls the voice leading in the case of ties between upperNote and lowerNote
		
		THIS STILL NEEDS TESTING!
		
		"""
		originalTarget = target 
		
		"""get the scale degree closest to our current value"""
		if recurse and self.container:
			target = self.container.scale_degree(target, direction, recurse = True)		# target is now a scale degree
		
		PC = target % abs(self.modulus)
		upperPC = min(self.wrappedSortedPCs, key = lambda x: x - PC if x >= PC else INFINITY)
		
		path = upperPC - PC
		upperScaleTone = target + path						
		
		sdIndex = self.PCs.index(upperPC % self.modulus)
		referenceTone = self.currentNotes[sdIndex]					# the right PC but could be in a different octave
		
		upperSD = self.SDs[sdIndex] + int((upperScaleTone - referenceTone)/self.modulus)*self.len		# calculate octave difference and add that number of degrees
		
		if path == 0 or (direction and direction > 0):
			return upperSD
			
		index = self.wrappedSortedPCs.index(upperPC)
		lowerScaleTone = self.sortedPCs[index - 1]					# NOT wrappedSortedPCs because the upperPC could be the first in the list!
		
		if lowerScaleTone > PC:
			lowerScaleTone = lowerScaleTone - self.modulus
		
		lowerScaleTone = target + lowerScaleTone - PC
		
		sdIndex2 = self.PCs.index(lowerScaleTone % self.modulus)
		referenceTone = self.currentNotes[sdIndex2]
		
		lowerSD = self.SDs[sdIndex2] + int((lowerScaleTone - referenceTone)/self.modulus)*self.len
		
		"""quantize down"""
		if direction and direction < 0:
			return lowerSD
		
		"""if we are recursing, we want to compare chromatic distances; otherwise scale-degree"""
		if recurse:
			target = originalTarget												# target is now a chromatic pitch
		
		frac1 = (target - self.__getitem__(lowerSD, recurse = recurse))			# (positive) distance to lower tone
		frac2 = (self.__getitem__(upperSD, recurse = recurse) - target)			# (positive) distance to upper tone
		
		self.quantizeSteps = [frac1, frac2]
		
		if direction is not None:
			"""quantize to the nearest scale tone"""
			if frac2 < frac1:
				return upperSD
			if frac1 < frac2:
				return lowerSD
			if tieBreak == None:								# randomly choose when confronted with ties
				return random.choice([upperSD, lowerSD])
			if tieBreak > 0:
				return upperSD
			return lowerSD
		
		"""
		if the nearest neighbors are adjacent scale degrees we can return a float as a scale degree;
			otherwise we need a list because a float makes no sense
			
		self.quantizeSteps contains the fractions
			
		"""
		
		if upperSD - lowerSD == 1:
			return lowerSD + (frac1 / (frac1 + frac2))
		return [lowerSD, upperSD]

class ScalePattern(Default):
	
	"""
	A ScalePattern groups, orders, and embellishes notes in a scale.  It is a component of a MultiPattern, which can manipulate multiple patterns.
	
	This is a rewriting of the Arpeggiation class, which requires a chord.  It is a first step toward eliminating chords, which are more trouble than they're worth.
	
	The main things it supplies are (a) ordering; and (b) NHT syntax.
	
	TODO: write "save" and "restore" routines that allow for cardinality changes.
	
	NHT syntax: 
	
		Incomplete neighbors (can also be used for passing tones or complete neighbors):
			['N', targetNote, hierarchicalLevel, interval]
				here the nonharmonic tone is targetNote + interval within the scale given by hierarchicalLevel
				hierarchicalLevel of 'S' represents a chordal skip
				TODO: control quantization direction, which should be the opposite of the skip
	
		Unanchored neighbors (displacements):
			['F', targetNote, targetLevel, interval, intervalLevel]
				basically a nonharmonic tone but not directed at another note in the motive
				targetLevel is relative to the chord's own scale, so one level deeper than most nonharmonic tones
	
		Inversional nonharmonic tones:
			['I', targetNote, hierarchicalLevel, indexNumber]

		Rests are nonharmonic tones, useful for arpeggiation:
			['R']

		Tacit harmonic notes (useful for initializing a motive that doesn't have every note in its chord)
			['T', 60]

		Pedal tones:
			['P', pedalNote, hierarchicalLevel, interval]
				pedalNote is a scale degree at hierarchical level1 
			TODO: these are unanchored neighbors, just translate them into that language
	
		TODO: transposable pedal, or addressable notes
			['CHORDELEMENT', 'R'] [root, set up a chordelement dictionary]
			
		TODO:
		Quantized value (Partian tintinabulation)
			['Q', targetNote, hierarchicalLevel, interval, quantizeDirection]
		
			['N', targetNote, hierarchicalLevel, interval, quantizationObject, quantizationDirection]
	
		Hierarchical levels are listed from surface to background:
	
			Level 0 ... x containing scales
			Level None chromatic scale
		
		These all need a text form as well:
			3N4–1 nonharmonic tone
	
	"""
	
	defaults = {'startNotes': [60, 62, 64],
							
							'modulus': 12,
							'multiset': False,
							'useMidi': True,
							
							'container': None,
							'chord': None,
							'outPCs': [],
							'newScale': False,
							
							'autoQuantize': {'direction': 0, 'tieBreak': 1},			# arguments to scale quantize
							
							'noteGroups': [],
							'fullStructure': [],
							'screenedNotes': [],
							'durations': [1],
							
							'objectDict': {},					# for multi-object structures
							
							'currentOutput': [],
							
							'useGroups': True,
							
							'chainMotive': False,
							'lastChordTone': False,
							
							'iterIndex': 0,
							'noteIndex': 0
						}
	
	def __init__(self, *args, **kwargs):
		
		self.set_default_arguments(kwargs)
			
		if args:
			self.startNotes = args[0]
			if type(self.startNotes) is str:
				self.startNotes = self.parse_string(self.startNotes)
			
		self.modulus = len(self.container) if self.container else 12
		
		self.initialStartNotes = self.startNotes[:]
		
		if self.startNotes:
			self.fullStructure, self.noteGroups, screenedNotes = self.parse_nonharmonic_tones()
			if not self.screenedNotes:
				self.screenedNotes = screenedNotes
			if not self.noteGroups:
				self.noteGroups = noteGroups
		
		if self.defaultArgs.get('useGroups'):				# should be irrelevant?
			self.useGroups = True
		
	def parse_string(self, s):
		"""
		Text syntax:
		
			brackets identify groups
		
			[atom1 atom2] []
		
			xANYTHING - tacit note
		
			xNy+z or xNy-z: neighbor tone at level x of the hierarchy, target note y, interval z (+ or - is required).  If x is missing it = 0.
			CSy+z or CSy-z: chordal skip, target y, steps z
			CNy+z or CNy-z: chromatic neighbor
		
			Px pedal tone x
		
		
		"""
		
		sSplit = parse_paren_unit(s.split(), '[', ']')
		sSplit = [x.replace('[', '').replace(']', '') for x in sSplit]
		
		sSplit = parse_paren_unit(sSplit, '(', ')')
		
		out = []
		totalSyms = 0
		lastNote = False
		
		for theSym in sSplit:
			symSplit = [x for x in theSym.split() if x]
			symSplit = parse_paren_unit(symSplit, '(', ')')
			if len(symSplit) == 1:
				targetList = out
				nestFlag = False
			else:
				targetList = []
				nestFlag = True
			for atom in symSplit:
				if atom.startswith('CS'):
					atom = atom.replace('CS', 'SN')
				quantizationData = []
				if atom.count('Q'):
					i = atom.index('Q')
					qData = atom[i+1:].split(' ')
					atom = atom[:i]
					varName = qData[0].replace('(', '')
					direction = int(parse_number(qData[1])) if len(qData) > 1 else 0
					tieBreak = int(parse_number(qData[2])) if len(qData) > 2 else 1
					extraTranspositions = int(parse_number(qData[3])) if len(qData) > 3 else 0
					quantizationData = [varName, direction, tieBreak, extraTranspositions]
				
				if atom.startswith('x') or atom.startswith('X'):
					tacitFlag = True
					atom = atom[1:]
				else:
					tacitFlag =  False
				if atom.count('N'):
					obj = self.parse_neighbor(atom) + quantizationData
					lastNote = False
				elif atom.startswith('P'):										# pedal tone, need to have level as a settable parameter
					obj = ['P', parse_lettername(atom[1:]), None]
					lastNote = copy.deepcopy(obj)
				elif atom == '-':												# shorthand for repeating the last note
					if lastNote:						
						obj = lastNote
					else:									
						obj = ['N', totalSyms - 1, None, 0]
				elif atom == 'R':
					obj = ['R']
					lastNote = copy.deepcopy(obj)
				else:
					obj = parse_lettername(atom)
					lastNote = obj
				if tacitFlag:
					obj = ['T', obj]
				targetList.append(obj)
				totalSyms += 1
			if nestFlag:
				out.append(targetList)
				lastNote = copy.deepcopy(targetList)
		return out
	
	def parse_neighbor(self, s):
		if s.count('('):
			i = s.index('(')
			j = s.index(')')
			objectString = s[i+1:j]
			s = s[:i] + s[j+1:]
			varName, num, *junk = objectString.split()
			target = (self.get_object(varName), int(num))
		else:
			target = None
		level, data = s.split('N')
		if data.count("+"):
			data = data.split('+')
			sign = 1
		else:
			data = data.split('-')
			sign = -1
		if not level:
			level = 0
		elif level.count('C'):
			level = None
		elif level.count('S'):
			level = 'S'
		else:
			level = int(level)
		if target is None:
			target = int(data[0])
		interval = sign*int(data[1])
		return ['N', target, level, interval]
		
	def get_object(self, varName):
		qObj = self.objectDict.get(varName)
		if not qObj:
			print("ARCA ERROR: Can't find quantization object:", varName)
			return False
		self.associatedObjects.add(qObj)
		return qObj
	
	def t(self, littleT = 0):
		self.scale.t(littleT)
		
	def quantize_MIDI_note(self, n):
		sd = self.container.scale_degree(n, direction = None, recurse = True)
		if type(sd) is int:
			return sd
		if not self.autoQuantize: 
			print("NONHARMONIC NOTE NOT IN CONTAINER, choosing upper neighbor", n)
			return int(round(sd, 0))
		return float(n)
	
	def autoquantize(self, n):
		n = int(n)
		for i, c in enumerate(self.scale.hierarchy):
			"""This is not tested"""
			testSD = c.scale_degree(n, direction = None, recurse = True)
			if type(testSD) is int:
				target = self.container.scale_degree(n, recurse = True, **self.autoQuantize)
				targetMIDI = self.container[target]
				targetAtLevel = c.scale_degree(targetMIDI, recurse = True)
				return ['F', target, 0, testSD - targetAtLevel, i]							# ['F', targetNote, targetLevel, interval, intervalLevel]
		target = self.scale.scale_degree(n, recurse = True, **self.autoQuantize)
		targetMIDI = self.container[target]
		return ['F', target, 0, n - targetMIDI, None]
		
	def parse_nonharmonic_tones(self, notes = None):
		
		if notes is None:
			notes = self.startNotes[:]
		
		notes, noteGroups = flatten_note_list(notes)
		totalNotes = []
		

		"""First, convert MIDI notes to scale degrees in the containing scale"""
		"""TODO: Test this to make sure it works with multisets, which should be handled at the scale level"""
		if self.useMidi and self.container:
			newNotes = []
			for n in notes:
				if type(n) is not list:
					newNotes.append(self.quantize_MIDI_note(n))
				else:
					if n[0] == 'T':
						n[1] = self.quantize_MIDI_note(n[1])			# tacit note
						newNotes.append(n)
					elif n[0] == 'F':					
						level = min(0, n[2] - 1)								
						n[1] = self.container.hierarchy[level].scale_degree(n[1], recurse = True)		# ['F', targetNote, targetLevel, interval, intervalLevel]
						if level == 0: newNotes.append(n)
					newNotes.append(n)
			notes = newNotes
		
		"""Now, collect harmonic notes so we can make the motive's intrinsic scale"""
		harmonicNotes = []
		screenedNotes = []										# notes not to play
		
		for i, n in enumerate(notes):
			if type(n) is list:
				if n[0] == 'T':
					screenedNotes.append(i)
					harmonicNotes.append(n[1])
				elif n[0] == 'F' and n[2] == 0:
					harmonicNotes.append(n[1])
			elif type(n) is int:
				harmonicNotes.append(n)
		
		if self.container:
			keyOffset = self.container.keyOffset
			harmonicNotes = [self.container.scalarPC_internal(x) for x in harmonicNotes]
		else:
			keyOffset = 60
			harmonicNotes = [x % self.modulus for x in harmonicNotes]
			
		harmonicNotes = sorted(set(harmonicNotes))
		
		self.scale = Scale(harmonicNotes, container = self.container, keyOffset = keyOffset, useMidi = False)
		
		fullStructure = notes[:]								# don't need to rename I don't think
		
		for i, n in enumerate(fullStructure):
			
			if type(n) is float:								# once we have the intrinsic scale, autoquantize the floats
				fullStructure[i] = self.autoquantize(n)
			elif type(n) is list:
				if n[0] == 'F' and n[2] == 0:
					n[1] = self.scale.scale_degree(n[1], recurse = False)
			else:
				fullStructure[i] = (self.scale.scale_degree(n, recurse = False),)
				
		self.lastFullStructure = copy.deepcopy(fullStructure)
		
		if screenedNotes:
			newNoteGroups = []
			for n in noteGroups:
				if type(n) is int:
					if n not in screenedNotes:
						newNoteGroups.append(n)
				elif type(n) is list:
					l = [x for x in n if x not in screenedNotes]
					if l:
						newNoteGroups.append(l)
			noteGroups = newNoteGroups
					
		return fullStructure, noteGroups, screenedNotes
			
	def evaluate_NHT(self, fullStructureIndex = 0, fullStructure = None): #note = 60, level = 0, alteration = lambda x: x):
		
		"""	
			This is a third version of evaluate_NHT that is nonrecursive and which builds chords note-by-note.  
			
			It can handle chord changes during the life of a motive.
		
			as you go through fullStructure starts as a specification of the motive, and gradually gets you
			which the recursive structures cannot.
		
			TODO: Incorporate this new logic into the CoreVoiceleading version of the routine.
		
			TODO: save all the scale degrees of the different NHTs?
		
			TODO: forward looking NHT evaluation?
		
		
		"""
		
		if fullStructure is None:
			fullStructure = self.fullStructure
		
		n = fullStructure[fullStructureIndex]
		
		"""tuples are scale degrees; integers are chromatic MIDI notes that have already been calculated"""
		if type(n) is tuple:
			return self.scale[n[0]]
			
		elif type(n) is int:
			return n
			
		elif n is None:
			return None
		
		NHTtype = n[0]
		
		if NHTtype == 'R':						# rest
			return None
			
		elif NHTtype == 'P':
			targetNote = n[1]
			level = n[2]
			"""TODO: test this, I suspect it doesn't work exactly"""
			if self.container and level is not None:
				targetNote = self.container.hierarchy[level].__getitem__(note)
			return self.scale.output_note(targetNote)
		
		elif NHTtype == 'F':					# floating nonharmonic tone ['F', targetNote, targetLevel, intervalLevel, interval]
			targetNote = self.scale.get_note_in_hierarchy(n[1], level = n[2])
			parameter = n[3]
			if n[4] != None:
				level = n[4] - 1
			else:
				level = None
			quantize = False
		
		else:
			"""first recurse through the NHT chain, and then the collectional chain, to get the chromatic representation of the target note"""
			if type(n[1]) is not int:
				targetnote = n[1][0].__getitem__(n[1][1])
			
			targetNote = self.evaluate_NHT(n[1], fullStructure)
			level = n[2]
			parameter = n[3]
			quantize = len(n) > 4 
		
		"""
		
		Test implementation of chordal skips
				TODO: allow quantization
		
		"""
		
		if level == 'S':
			""" TODO: allow control of quantization """
			
			if len(n) > 4:
				quantizeDirection = n[4]
			else:
				quantizeDirection = -parameter
				
			targetNote = self.scale.quantize(targetNote, direction = quantizeDirection)
			
			return self.scale.output_note(targetNote + parameter)
			
		"""back up the hierarchy to the NHT level"""
		if self.container:
			if level and level >= len(self.container.hierarchy):
				level = None
			if level is not None:
				level = level % len(self.container.hierarchy)
				targetNote = self.container.hierarchy[level].scale_degree(targetNote, direction = 0)			# find the scale degree at the appropriate level corresponding to the chromatic note
		else:
			level = None
		
		"""make the alteration"""
		if NHTtype == 'I':
			targetNote = parameter - targetNote
		else:
			targetNote = targetNote + parameter
		
		"""back down the collectional hierarchy to the chromatic scale"""
		if level is not None:
			targetNote = self.container.hierarchy[level].__getitem__(targetNote, recurse = True)
		
		"""if we want to quantize, do it here"""
		if quantize:								
			qObj = self.objectDict.get(n[4])
			direction = n[5] if len(n) > 5 else 0
			tieBreak = n[6] if len(n) > 6 else 1
			extraTransposition = n[7] if len(n) > 7 else 0
			targetNote = qObj.quantize(targetNote, direction = direction, tieBreak = tieBreak, extraTransposition = extraTransposition)	
		
		return self.scale.output_note(targetNote)
	
	def get_durations(self):
		durations = []
		for i in range(len(self.iterObject)):
			durations.append(self.durations[(self.noteIndex + i)%len(self.durations)])
		self.noteIndex += i + 1
		return durations
		
	def __iter__(self):
		self.iterIndex = -1
		return self
	
	def __next__(self):
		self.iterIndex += 1
		if self.iterIndex < len(self.iterObject):
			if self.useGroups:
				return self.get_group(self.iterIndex)
			return self[self.iterIndex]
		else:
			raise StopIteration
		
	def __len__(self):
		return len(self.iterObject)
	
	def __getitem__(self, key): #, makeList = True):
		
		if isinstance(key, slice):
			start = key.start if (key.start is not None) else 0
			end = key.stop if (key.stop is not None) else len(self.fullStructure)
			step = key.step if (key.step is not None) else 1
			return [self.build_output(x) for x in range(start, end, step)]
		elif isinstance(key, list):
			return [self.build_output(x) for x in key]
		
		return self.build_output(key)
	
	def invert_motive(self):
		for l in self.fullStructure:
			if l[0] in 'DNF':
				l[3] = l[3]*-1
	
	def transform(self, s):
		s = s.replace(' ', '')
		if s.startswith('*'):
			return self.scale.transform(s[1:])
		elif s.startswith('i'):
			self.invert_motive()
		else:
			return self.scale.transform(s)
	
	def get_notes(self):
		return [self.evaluate_NHT(x) for x in [i for i in range(len(self.fullStructure)) if i not in self.screenedNotes]]
	
	@property
	def useGroups(self):
		if self.iterObject == self.noteGroups:
			return True
		return False
	
	@useGroups.setter
	def useGroups(self, val):
		if val:
			self.iterObject = self.noteGroups
		else:
			self.iterObject = [x for x in range(len(self.fullStructure)) if x not in self.screenedNotes]
	
	def build_output(self, i = 0, reset = False):
		
		if reset or i == 0 or (not self.currentOutput):
			self.currentOutput = self.fullStructure[:]
			if self.chainMotive and self.lastChordTone:
				self.currentOutput[0] = copy.deepcopy(self.lastChordTone)
		
		if i != 0 and (not (type(self.currentOutput[i-1]) is int)):
			self.build_output(i-1)
		
		n = self.evaluate_NHT(i, self.currentOutput)
		self.currentOutput[i] = n
		
		"""n is chromatic, sd is a scale degree"""
		"""TODO: fix this, """
		if n and self.chainMotive:
			self.lastChordTone = n
			print("ERROR: CHAINMOTIVE NEEDS TO BE IMPLEMENTED FOR SCALEPATTERN")
			"""if self.container:
				sd = self.container.scale_degree(n)
			ct = self.chord.get_chord_tone(sd)
			if ct:
				self.lastChordTone = ct"""
				
		return n
	
	def build_groups(self, groupNumber = 0, reset = None):
		if reset is None:
			reset = (groupNumber == 0)
		index = self.noteGroups[groupNumber]
		if type(index) is int:
			self.build_output(index, reset = reset)
			return [self.currentOutput[index]]
		if type(index) is list:
			self.build_output(max(index), reset = reset)
			return [self.currentOutput[x] for x in index]
		self.build_output(index.stop - 1, reset = reset)
		return self.currentOutput[index]
	
	def get_group(self, index):
		return self.build_groups(index)
	
	@property
	def notes(self):
		return self.get_notes()

class Arpeggiation(Default):
	
	"""

	
	An arpeggiation pattern groups, orders, and embellishes the notes of a chord.
	
	The main thing it supplies is NHT syntax [though this should probably be moved to CoreVoiceLeading]
	
	NHT syntax: 
	
		Incomplete neighbors (can also be used for passing tones or complete neighbors):
			['N', targetNote, hierarchicalLevel, interval]
				here the nonharmonic tone is targetNote + interval within the scale given by hierarchicalLevel
				hierarchicalLevel of 'S' represents a chordal skip
	
		Unanchored neighbors:
			['F', targetNote, targetLevel, interval, intervalLevel]
				basically a nonharmonic tone but not directed at another note in the motive
				targetLevel is relative to the chord's own scale
	
		Pedal tones:
			['P', pedalNote, hierarchicalLevel, interval]
				pedalNote is a scale degree at hierarchical level1 
				TODO: these are bascially unanchord neighbors, they can be made obsolete
		
		TODO: These are obsolete
		Octave doublings: ['D', targetNote, hierarchicalLevel, interval]
			neighbor tones created by Chord, representing note Duplications.  
			Functionally equivalent to neighbors, but updated differently.
	
		TODO: transposable pedal, or addressable notes
			['CHORDELEMENT', 'R'] [root, set up a chordelement dictionary]

		Or for inversional nonharmonic tones:
			['I', targetNote, hierarchicalLevel, indexNumber]

		Rests are nonharmonic tones, useful for arpeggiation:
			['R']

		Tacit harmonic notes (useful for initializing a motive that doesn't have every note in its chord)
			['T', 60]
			
		TODO:
	
		Quantized value (Partian tintinabulation)
			['Q', targetNote, hierarchicalLevel, interval, quantizeDirection]
		
			['N', targetNote, hierarchicalLevel, interval, quantizationObject, quantizationDirection]
	
		For most NHTs, herarchical levels are listed from surface to background:
	
			Level 0 ... x containing scales
			Level None chromatic scale
	
		For the Unanchored Neighbors, the scale of the motive is level 0, so all those numbers are increased by one.
		
		These all need a text form as well:
			3N4–1 nonharmonic tone
	
	"""
	
	defaults = {'startNotes': [60, 62, 64],
							
							'modulus': 12,
							'multiset': False,
							'useMidi': True,
							
							'container': None,
							'outPCs': [],
							'newChord': False,
							
							'currentNotes': [],
							'noteGroups': [],
							'fullStructure': [],
							'screenedNotes': [],
							'durations': [1],
							
							'objectDict': {},					# for multi-object structures
							'dependentObjects': [],
							
							'currentOutput': [],
							
							'useGroups': True,
							
							'chainMotive': False,
							'lastChordTone': False,
							
							'iterIndex': 0,
							'noteIndex': 0
							
						}
	
	def __init__(self, *args, **kwargs):
		
		self.set_default_arguments(kwargs)
			
		if args:
			self.startNotes = args[0]
			if type(self.startNotes) is str:
				self.startNotes = self.parse_string(self.startNotes)
			if len(args) > 1:
				self.chord = args[1]
			
		if self.chord is not None:
			self.container = self.chord.container
			self.modulus = self.chord.modulus
			self.multiset = self.chord.multiset
		else:
			self.modulus = len(self.container) if self.container else 12
		
		self.initialStartNotes = self.startNotes[:]
		
		if self.startNotes:
			if self.newChord or self.chord is None:
				self.currentNotes, self.fullStructure, self.noteGroups, self.screenedNotes = self.parse_nonharmonic_tones()
			else:
				junk, self.fullStructure, self.noteGroups, self.screenedNotes = self.parse_nonharmonic_tones_without_replacement()
				self.currentNotes = self.chord.currentNotes
		
		if self.chord is None:
			useMidi = False if self.container else self.container
			self.chord = Chord(self.currentNotes[:], mainPattern = self, container = self.container, multiset = self.multiset, useMidi = useMidi)
		
		self.chord.dependentObjects.append(self)
			
		if self.defaultArgs.get('useGroups'):
			self.useGroups = True
		
	def parse_string(self, s):
		"""
		Text syntax:
		
			brackets identify groups
		
			[atom1 atom2] []
		
			xANYTHING - tacit note
		
			xNy+z or xNy-z: neighbor tone at level x of the hierarchy, target note y, interval z (+ or - is required).  If x is missing it = 0.
			CSy+z or CSy-z: chordal skip, target y, steps z
			CNy+z or CNy-z: chromatic neighbor
		
			Px pedal tone x
		
		
		"""
		
		sSplit = parse_paren_unit(s.split(), '[', ']')
		sSplit = [x.replace('[', '').replace(']', '') for x in sSplit]
		
		sSplit = parse_paren_unit(sSplit, '(', ')')
		
		out = []
		totalSyms = 0
		lastNote = False
		
		for theSym in sSplit:
			symSplit = [x for x in theSym.split() if x]
			symSplit = parse_paren_unit(symSplit, '(', ')')
			if len(symSplit) == 1:
				targetList = out
				nestFlag = False
			else:
				targetList = []
				nestFlag = True
			for atom in symSplit:
				if atom.startswith('CS'):
					atom = atom.replace('CS', 'SN')
				quantizationData = []
				if atom.count('Q'):
					i = atom.index('Q')
					qData = atom[i+1:].split(' ')
					atom = atom[:i]
					varName = qData[0].replace('(', '')
					direction = int(parse_number(qData[1])) if len(qData) > 1 else 0
					tieBreak = int(parse_number(qData[2])) if len(qData) > 2 else 1
					extraTranspositions = int(parse_number(qData[3])) if len(qData) > 3 else 0
					quantizationData = [varName, direction, tieBreak, extraTranspositions]
				
				if atom.startswith('x') or atom.startswith('X'):
					tacitFlag = True
					atom = atom[1:]
				else:
					tacitFlag =  False
				if atom.count('N'):
					obj = self.parse_neighbor(atom) + quantizationData
					lastNote = False
				elif atom.startswith('P'):										# pedal tone, need to have level as a settable parameter
					obj = ['P', parse_lettername(atom[1:]), None]
					lastNote = copy.deepcopy(obj)
				elif atom == '-':												# shorthand for repeating the last note
					if lastNote:						
						obj = lastNote
					else:									
						obj = ['N', totalSyms - 1, None, 0]
				elif atom == 'R':
					obj = ['R']
					lastNote = copy.deepcopy(obj)
				else:
					obj = parse_lettername(atom)
					lastNote = obj
				if tacitFlag:
					obj = ['T', obj]
				targetList.append(obj)
				totalSyms += 1
			if nestFlag:
				out.append(targetList)
				lastNote = copy.deepcopy(targetList)
		return out
	
	def parse_neighbor(self, s):
		if s.count('('):
			i = s.index('(')
			j = s.index(')')
			objectString = s[i+1:j]
			s = s[:i] + s[j+1:]
			varName, num, *junk = objectString.split()
			target = (self.get_object(varName), int(num))
		else:
			target = None
		level, data = s.split('N')
		if data.count("+"):
			data = data.split('+')
			sign = 1
		else:
			data = data.split('-')
			sign = -1
		if not level:
			level = 0
		elif level.count('C'):
			level = None
		elif level.count('S'):
			level = 'S'
		else:
			level = int(level)
		if target is None:
			target = int(data[0])
		interval = sign*int(data[1])
		return ['N', target, level, interval]
	
	def get_object(self, varName):
		qObj = self.objectDict.get(varName)
		if not qObj:
			print("ARCA ERROR: Can't find quantization object:", varName)
			return False
		self.associatedObjects.add(qObj)
		return qObj
	
	def transpose_along_chord(self, noteList, littleT = 0):
		tDict = self.chord.make_Tt_dict(0, littleT)
		"""TODO: need more here"""
	
	def t(self, littleT = 0):
		
		tDict = self.chord.make_Tt_dict(0, littleT)
		myNotes = copy.deepcopy(self.fullStructure)
		
		for i, n in enumerate(myNotes):
			if type(n) is list:
				if n[0] == 'D':
					targetNote = self.chord.currentNotes[n[1]]
					interval = int(self.modulus*n[3]/12)
					myNotes[i] = targetNote + interval + tDict.get(targetNote % self.modulus, 0)
		
		output = self.parse_nonharmonic_tones_without_replacement(myNotes, useMidi = False)
		junk, self.fullStructure, junk, junk = output
		
	def parse_nonharmonic_tones(self, notes = None):
		
		if notes is None:
			notes = self.startNotes[:]
		
		notes, noteGroups = flatten_note_list(notes)
		
		if self.useMidi and self.container:
			newNotes = []
			for n in notes:
				if type(n) is not list:
					newNotes.append(self.container.scale_degree(n, recurse = True))
				else:
					if n[0] == 'T' and type(n[1]) is not list:
						newNotes.append(['T', self.container.scale_degree(n[1], recurse = True)])
					else:
						newNotes.append(n)
			notes = newNotes
		
		
		strippedNotes = []										# only the harmonic tones
		strippedPCs = []										# the same but their PCs
		pcCounts = {}											# the count of the PCs
		fullStructure = notes[:]								# the input, but with notenumbers replaced by indexes
		screenedNotes = []										# notes not to play
		
		for i, n in enumerate(notes):
			
			if type(n) is list and len(n) >= 2 and n[0] == 'T':
				n = n[1]
				fullStructure[i] = n
				screenedNotes.append(i)
				
			if type(n) is not list:
				PC = n % self.modulus
				if type(PC) is float:
					PC = round(PC, 2)
					n = round(n, 2)
				count = pcCounts.get(PC, -1) + 1
				pcCounts[PC] = count
				
				"""
				
				There are two ways to treat duplicate PCs:
					as nonharmonic notes [octave doublings]
					or as multiests (fractionally displaced)
				
				"""
				if count:
					if not self.multiset:									# convert duplicate to nonharmonic tone
						index = strippedPCs.index(PC)
						item = strippedNotes[index]
						fullStructure[i] = ['D', strippedNotes.index(item), None, int(12*((n-item)/self.modulus))]
						continue
					else:
						n += .0001*count									# turn duplicate into fractionally displaced PC
				fullStructure[i] = ['D', len(strippedNotes), None, 0]
				strippedNotes.append(n)
				strippedPCs.append(PC)
		
		self.lastFullStructure = copy.deepcopy(fullStructure)
		
		"""
		if self.makeRecursiveNeighbors:
			for i, n in enumerate(fullStructure):
				fullStructure[i] = self.create_recursive_neighbors(n, fullStructure)"""
		
		if screenedNotes:
			newNoteGroups = []
			for n in noteGroups:
				if type(n) is int:
					if n not in screenedNotes:
						newNoteGroups.append(n)
				elif type(n) is list:
					l = [x for x in n if x not in screenedNotes]
					if l:
						newNoteGroups.append(l)
			noteGroups = newNoteGroups
					
		return strippedNotes, fullStructure, noteGroups, screenedNotes
			
	def parse_nonharmonic_tones_without_replacement(self, notes = None, useMidi = None):			# problem passing None to brython?
		
		
		"""
		PARSE NONHARMONIC TONES BUT PRESERVE PCS
			Used when multiple ArpeggiationPatterns use the same chord
			
		"""
		
		if notes is None:
			notes = self.startNotes
		
		notes, noteGroups = flatten_note_list(notes)
		
		if useMidi or (useMidi is None and self.useMidi and self.container):
			newNotes = []
			for n in notes:
				if type(n) is not list:
					newNotes.append(self.container.scale_degree(n, recurse = True))
				else:
					if n[0] == 'T' and type(n[1]) is not list:
						newNotes.append(['T', self.container.scale_degree(n[1], recurse = True)])
					else:
						newNotes.append(n)
			notes = newNotes
		
		strippedNotes = self.chord.currentNotes											# only the harmonic tones
		strippedPCs = self.chord.PCs													# have to use the order of strippedNotes; self.chord.sortedPCs is sorted
		fullStructure = notes[:]														# the input, but with notenumbers replaced by indexes
		screenedNotes = []																# notes not to play
		
		for i, n in enumerate(notes):
			
			if type(n) is list and len(n) >= 2 and n[0] == 'T':
				n = n[1]
				fullStructure[i] = n
				screenedNotes.append(i)
				
			if type(n) is not list:
				PC = n % self.modulus
				if type(PC) is float:
					if self.multiset:
						PC = round(PC, 4)
					else:
						PC = round(PC, 2)
					n = round(n, 2)
				
				if PC not in strippedPCs:
					print("ARCA ERROR: can't find chord tone for note", n)
					print(strippedNotes, strippedPCs)
					print(self.chord.notes)
					continue
				
				index = strippedPCs.index(PC)
				item = strippedNotes[index]
				
				interval = int(12*((n-item)/self.modulus))
				
				fullStructure[i] = ['D', index, None, interval]
		
		self.lastFullStructure = copy.deepcopy(fullStructure)
		
		if self.makeRecursiveNeighbors:
			for i, n in enumerate(fullStructure):
				fullStructure[i] = self.create_recursive_neighbors(n, fullStructure, updateTargets = False)
		
		if screenedNotes:
			newNoteGroups = []
			for n in noteGroups:
				if type(n) is int:
					if n not in screenedNotes:
						newNoteGroups.append(n)
				elif type(n) is list:
					l = [x for x in n if x not in screenedNotes]
					if l:
						newNoteGroups.append(l)
			noteGroups = newNoteGroups
		
		return strippedNotes, fullStructure, noteGroups, screenedNotes
	
	def create_recursive_neighbors(self, n, fullStructure, updateTargets = True):
		
		"""
		TODO: obsolete!
		
		We might want neighbor notes to be self-standing and not involve any reference to fullStructure for its evaluation; 
		this is because we might want neighbors that are self-standing, in the case where the Chord is embedded in a larger object such as a multicollection PitchMotive
		
		This turns out not to work very well, because chords can change during the life of a motive
		
		"""
		
		if type(n) is list and len(n) >= 4:
			if type(n[1]) is list: return n
			target = fullStructure[n[1]]
			if type(target) is not list:
				if updateTargets or n[0] == 'N':
					return [n[0], target] + n[2:]														# replace fullStructure index with currentNotes index
				else:
					return n[:]
			else:
				if n[0] == 'D':
					return n[:]
				if target[0] == 'D':
					return [n[0], target] + n[2:]	
				return [n[0], self.create_recursive_neighbors(target, fullStructure)] + n[2:]		# create recursive NHTs
		return n
			
	def evaluate_NHT(self, fullStructureIndex = 0, fullStructure = None): #note = 60, level = 0, alteration = lambda x: x):
		
		"""	
			This is a third version of evaluate_NHT that is nonrecursive and which builds chords note-by-note.  
			
			It can handle chord changes during the life of a motive.
		
			as you go through fullStructure starts as a specification of the motive, and gradually gets you
			which the recursive structures cannot.
		
			TODO: Incorporate this new logic into the CoreVoiceleading version of the routine.
		
			TODO: save all the scale degrees of the different NHTs?
		
			TODO: forward looking NHT evaluation?
		
		I fear the logic here is still not geat.
		
		"""
		
		if fullStructure is None:
			fullStructure = self.fullStructure
		
		n = fullStructure[fullStructureIndex]
		
		"""ints are chromatic MIDI notes that have already been calculated"""
		if type(n) is int:
			return n
		
		NHTtype = n[0]
		
		if NHTtype == 'R':						# rest
			return None
			
		elif NHTtype == 'P':
			targetNote = n[1]
			level = n[2]
			"""TODO: test this, I suspect it doesn't work exactly"""
			if self.container and level is not None:
				targetNote = self.container.hierarchy[level].__getitem__(note)
			return self.chord.output_note(targetNote)
			
		elif NHTtype == 'D':					# octave doubling, should be changed to 'O'?
			
			targetNote = self.chord[n[1]]
			if type(targetNote) is list:
				targetNote = targetNote[0]
			return self.chord.output_note(targetNote + n[3])
		
		elif NHTtype == 'F':					# floating nonharmonic tone ['F', targetNote, targetLevel, interval, intervalLevel]
			targetNote = self.chord.scale.get_note_in_hierarchy(n[1], level = n[2])
			parameter = n[3]
			if n[4] != None:
				level = n[4] - 1
			else:
				level = None
			quantize = False
		
		else:
			"""first recurse through the NHT chain, and then the collectional chain, to get the chromatic representation of the target note"""
		
			if type(n[1]) is not int:
				targetnote = n[1][0].__getitem__(n[1][1])
			
			targetNote = self.evaluate_NHT(n[1], fullStructure)
			level = n[2]
			parameter = n[3]
			quantize = len(n) > 4 
		
		"""
		
		Test implementation of chordal skips
				TODO: allow quantization
		
		"""
		
		if level == 'S':
			""" TODO: allow control of quantization """
			
			if len(n) > 4:
				quantizeDirection = n[4]
			else:
				quantizeDirection = -parameter
				
			targetNote = self.chord.quantize(targetNote, direction = quantizeDirection)
			
			if self.container:
				targetNote = self.container.scale_degree(targetNote)
			
			tDict = self.chord.make_Tt_dict(0, parameter)
			targetNote = targetNote + tDict.get(targetNote % self.modulus, 0)
			
			if self.container:
				targetNote = self.container.__getitem__(targetNote, recurse = True)
			
			return self.chord.output_note(targetNote)
			
		"""back up the hierarchy to the NHT level"""
		if self.container:
			if level and level >= len(self.container.hierarchy):
				level = None
			if level is not None:
				level = level % len(self.container.hierarchy)
				targetNote = self.container.hierarchy[level].scale_degree(targetNote, direction = 0)			# find the scale degree at the appropriate level corresponding to the chromatic note
		else:
			level = None
		
		"""make the alteration"""
		if NHTtype == 'I':
			targetNote = parameter - targetNote
		else:
			targetNote = targetNote + parameter
		
		"""back down the collectional hierarchy to the chromatic scale"""
		if level is not None:
			targetNote = self.container.hierarchy[level].__getitem__(targetNote, recurse = True)
		
		"""if we want to quantize, do it here"""
		if quantize:								
			qObj = self.objectDict.get(n[4])
			direction = n[5] if len(n) > 5 else 0
			tieBreak = n[6] if len(n) > 6 else 1
			extraTransposition = n[7] if len(n) > 7 else 0
			targetNote = qObj.quantize(targetNote, direction = direction, tieBreak = tieBreak, extraTransposition = extraTransposition)	
		
		return self.chord.output_note(targetNote)
	
	def get_durations(self):
		durations = []
		for i in range(len(self.iterObject)):
			durations.append(self.durations[(self.noteIndex + i)%len(self.durations)])
		self.noteIndex += i + 1
		return durations
		
	def __iter__(self):
		self.iterIndex = -1
		return self
	
	def __next__(self):
		self.iterIndex += 1
		if self.iterIndex < len(self.iterObject):
			if self.useGroups:
				return self.get_group(self.iterIndex)
			return self[self.iterIndex]
		else:
			raise StopIteration
		
	def __len__(self):
		return len(self.iterObject)
	
	def __getitem__(self, key): #, makeList = True):
		
		if isinstance(key, slice):
			start = key.start if (key.start is not None) else 0
			end = key.stop if (key.stop is not None) else len(self.fullStructure)
			step = key.step if (key.step is not None) else 1
			return [self.build_output(x) for x in range(start, end, step)]
		elif isinstance(key, list):
			return [self.build_output(x) for x in key]
		
		return self.build_output(key)
	
	def invert_motive(self):
		for l in self.fullStructure:
			if l[0] in 'DN':
				l[3] = l[3]*-1
	
	def transform(self, s):
		s = s.replace(' ', '')
		if s.startswith('*'):
			return self.chord.transform(s[1:])
		elif s.startswith('t'):
			s = ''.join([x for x in s if x.isdigit() or x == '-'])
			i = int(s)
			self.t(i)
		elif s.startswith('i'):
			self.invert_motive()
	
	def get_notes(self):
		return [self.evaluate_NHT(x) for x in [i for i in range(len(self.fullStructure)) if i not in self.screenedNotes]]
	
	@property
	def useGroups(self):
		if self.iterObject == self.noteGroups:
			return True
		return False
	
	@useGroups.setter
	def useGroups(self, val):
		if val:
			self.iterObject = self.noteGroups
		else:
			self.iterObject = [x for x in range(len(self.fullStructure)) if x not in self.screenedNotes]
	
	def build_output(self, i = 0, reset = False):
		
		if reset or i == 0 or (not self.currentOutput):
			self.currentOutput = self.fullStructure[:]
			if self.chainMotive and self.lastChordTone:
				self.currentOutput[0] = copy.deepcopy(self.lastChordTone)
		
		if i != 0 and (not (type(self.currentOutput[i-1]) is int)):
			self.build_output(i-1)
		
		n = self.evaluate_NHT(i, self.currentOutput)
		self.currentOutput[i] = n
		
		"""n is chromatic, sd is a scale degree, ct is a chord tone of the form ['D', currentNotesIndex, 0, semitonalinterval]"""
		if n and self.chainMotive:
			if self.container:
				sd = self.container.scale_degree(n)
			ct = self.chord.get_chord_tone(sd)
			if ct:
				self.lastChordTone = ct
				
		return n
	
	def build_groups(self, groupNumber = 0, reset = None):
		if reset is None:
			reset = (groupNumber == 0)
		index = self.noteGroups[groupNumber]
		if type(index) is int:
			self.build_output(index, reset = reset)
			return [self.currentOutput[index]]
		if type(index) is list:
			self.build_output(max(index), reset = reset)
			return [self.currentOutput[x] for x in index]
		self.build_output(index.stop - 1, reset = reset)
		return self.currentOutput[index]
	
	def get_group(self, index):
		return self.build_groups(index)
	
	@property
	def notes(self):
		return self.get_notes()
		
class Chord(CoreVoiceleading):
	
	"""
	
	WARNING: DEPRECATED IN FAVOR OF SCALEPATTERN
	
	A Chord is a collection of harmonic notes with no pitch-class duplicatons (currentNotes) 
	plus one or more Arpeggiations (mainPattern) which can have PC doublings and nonharmonic tones.
	
	The CoreVoiceleading operations work on chords without PC duplications and without nonharmonic tones.
	
	Chord provides two ways of dealing with doublings:
	
		1. turning them into nonharmonic tones 
		2. creating multisets by displacing them by inaudible fractions (multset = tTrue)
	
	Arpeggiations are what I call "structured arpeggiations"
	
	TODO:
			we want a .sort function for the chord, so we need to register the fullStructures for updating
	
	TODO: 
		be better about passing arguments to mainPattern.  These should go by default.
	
	can initialize with either chordList = [...] for a round/repeating progression, or sequence = [...] for a sequence.  These are not default args
	
	
	"""
	defaults = {	'startNotes': [60, 62, 64],
							
							'modulus': 12,
							'container': None,
							'useMidi': True,
							'useGroups': True,
							'multiset': False,
							'outPCs': [],
							'mainPattern': False,
							'currentNotes': [],
							'objectDict': {},
							
							'durations': [1],
							
							'noteScreen': [],
							'dependentObjects': []
											
						}
	
	def __init__(self, *args, **kwargs):
		
		self.initialize_core()
		
		self.set_default_arguments(kwargs)
		
		if args:
			self.startNotes = args[0]
		
		self.initialize()
		
	def initialize(self):
		
		"""
		there are two ways to initialize a Chord: 
			using the container as its scale
			having it generate its own scale
		
		"""
		if self.container:
			if type(self.container) is int:
				self.container = Scale(self.container)
			self.modulus = self.container.len
			self.hierarchy = self.container.get_hierarchy()
			
		else:	
			self.modulus = 12
			self.hierarchy = []
		
		if self.multiset:
			self.output_note = self.output_note_multiset
		else:
			self.output_note = self.output_note_regular
			
		if hasattr(self, 'chordList') or hasattr(self, 'sequence'):
			if hasattr(self, 'sequence'):
				self.startNotes, prog = self.chords_to_transformations(self.sequence, schema = True)
			else:
				self.startNotes, prog = self.chords_to_transformations(self.chordList, schema = False)
			self.progression = prog
		
		
		if not self.mainPattern:
			self.mainPattern = Arpeggiation(self.startNotes, self, multiset = self.multiset, useMidi = self.useMidi, durations = self.durations, modulus = self.modulus, useGroups = self.useGroups, newChord = True, objectDict = self.objectDict)
			#self.durations = self.mainPattern.durations				# does this casuse error
		
		self.currentNotes = self.mainPattern.currentNotes[:]
		self.lastNotes = self.currentNotes[:]
		self.initialNotes = self.currentNotes[:]
		
		self.update()
		self.get_basic_VL()
	
	def update(self):
		self.PCs = [x % self.modulus for x in self.currentNotes]
		self.sortedPCs = sorted(self.PCs)
		self.make_scale()
	
	def get_chord_tone(self, n):
		PC = n % self.modulus
		try:
			index = self.PCs.index(PC)
		except ValueError:
			return None
		item = self.currentNotes[index]
		return ['D', index, None, int(12*((n-item)/self.modulus))]
		
	def quantize(self, *args, **kwargs):
		
		"""
		Quantizes to the current chord
		
		"""
		return self.scale.quantize(*args, **kwargs)
	
	def recurse_through_NHT_lists(self, l):
		if type(l) is list and len(l) >= 4 and type(l[1]) is list:
			return self.recurse_through_NHT_lists(l[1])
		return l
		
	def get_note_by_index(self, i):
		
		return self.mainPattern[i]
	
	def get_group(self, index):
		return self.mainPattern.get_group(index)
	
	def get_durations(self):
		return self.mainPattern.get_durations()
	
	def get_notes(self):
		return self.mainPattern.get_notes()
	
	def spell(self):
		return ' '.join(spell_MIDI(self.mainPattern.get_notes(), octaves = True))
		
	def apply_transposition_and_permutation(self, permDict = {}, transformDict = {}, targetList = None, updateNeighbors = True, updatePedals = False, invert = False):
		
		"""
		
		permutes and transforms targetList, with options to invert instead of transpose, affect NHTs, and apply to transposable pedals
		
		transformDict is either a collection of transpositions or a collection of index numbers 
		
		"""
		
		sign = -1 if invert else 1
		targetList = self.fullStructure if targetList is None else targetList	 
		out = [None]*len(targetList)
		
		for i, item in enumerate(targetList):
			targetIndex = permDict.get(i, i)
			targetTransp = transformDict.get(i, 0)
			if type(item) is not list:
				"""if item is not a list, it is an index to self.currentNotes"""
				self.currentNotes[item] = targetTransp + sign*self.currentNotes[item]
			else:
				"""item could be a recursively nested series of lists"""
				targetList = self.recurse_through_NHT_lists(item)
				if updatePedals and item[0] == 'Q':
					targetList[1] = targetTransp + sign*targetList[1]
				if updateNeighbors and item[0] == 'N':
					targetList[1] = targetIndex
					targetList[3] *= sign
			out[targetIndex] = item
				
		return out
	
	def get_permuted_item(self, item, originalVoiceNumber, permDict = {}, updateNeighbors = True):
		
		"""
		permutes and transforms targetList, with options to invert instead of transpose, affect NHTs, and apply to transposable pedals
		
		transformDict is either a collection of transpositions or a collection of index numbers
		
		"""
		
		targetIndex = permDict.get(originalVoiceNumber, originalVoiceNumber)
		
		if type(item) is not list:
			"""if item is not a list, it is an index to self.currentNotes"""
			return targetIndex
		else:
			"""item could be a recursively nested series of lists"""
			targetList = self.recurse_through_NHT_lists(item)
			if updateNeighbors and item[0] == 'N':
				targetList[1] = targetIndex
				targetList[3] *= sign
		return item
	
	def make_scale(self):
		self.scale = Scale(sorted(self.raw_content()), container = self.container, useMidi = True)
		return self.scale
	
	"""TODO: add recurse and depth parameters here?"""
	def __getitem__(self, key):
		
		if isinstance(key, slice):
			start = key.start if (key.start is not None) else 0
			end = key.stop if (key.stop is not None) else start + len(self.currentNotes)
			step = key.step if (key.step is not None) else 1
			return [self.__getitem__(x) for x in range(start, end, step)]
		
		if self.container:
			return self.container[self.currentNotes[key]]
		else:
			return self.currentNotes[key]
	
	def __contains__(self, item):
		return item in self.notes
		
	def MIDI(self):
		return self.get_notes()
	
	@property
	def core(self):
		return self.currentNotes
		
	@property
	def content(self):
		return self.get_notes()
	
	@property
	def notes(self):
		return self.get_notes()

class Pattern(Default):
	
	"""
	A Pattern is an arpeggiation structure containing a collection of chords, along with routines for keeping them consistent
		Pass it an arpeggiation pattern of the form [(Object, key), ..., (ChordObject, key)], withs sublists containing groups
		Key can be either an integer or slice, and will be passed to the Object's __getitem__method
		No Object will produce chromatic notes
	
		.sequence will cause the other chords to update rigidly
			- NOT with the same Tt, but with the correct transposition (or inversion)
			- referenceObject, if not specified, is the object with the largest currentNotes
	
		.broadcast will apply the same transformation to every chord (Tt or whatever)
	
	TODO: this needs a lot of work!!!
	
	"""
	
	defaults = {		'chordObjects': set(),
							
							'modulus': 12,
							'container': None,
							
							'useMidi': True,
							
							'mainObject': None,
							'array': []
											
						}
	
	def __init__(self, *args, **kwargs):
		
		self.set_default_arguments(kwargs)
		
		if args:
			self.startNotes = args[0]
			self.fullStructure, self.noteGroups, self.separatedPatterns = self.process_pattern(self.startNotes)
			self.chordObjects = set(self.separatedPatterns)
	
	@staticmethod
	def ignore_function(item):
		if type(item[0]) is not list:
			return True
		return False
	
	def process_pattern(self, inputNotes):
		flatList, noteGroups = flatten_note_list(inputNotes, ignoreFunction = Pattern.ignore_function)
		separatedPatterns = {}
		
		for chordObject, n in flatList:
			separatedPatterns.setdefault(chordObject, []).append(n)
		
		return flatList, noteGroups, separatedPatterns
	
	def next(self):
		for o in self.chordObjects:
			if o:
				o.next()
	
	def get_notes(self):
		out = []
		for o, n in self.fullStructure:
			if o:
				temp = o[n]
				if type(temp) is list:
					out += temp
				else:
					out.append(temp)
			else:
				out.append(n)
		return out
	
	def make_array(self):
		
		if self.array and self.mainObject:
			return
			
		self.mainObject = max(self.chordObjects, key = lambda x: len(x.currentNotes))
		theLen  = len(self.mainObject.currentNotes)
		self.array = []
		
		for c in self.chordObjects:
			
			if len(c.currentNotes) == theLen:
				self.array.append([c, list(range(len(c.currentNotes)))])
				
			elif len(c.fullStructure) == theLen:
				self.array.append([c, copy.deepcopy(c.fullStructure)])
				
			else:
				diff = theLen - len(c.currentNotes)
				self.array.append([c, list(range(theLen)) + ['R'] * diff])
					
		"""
		l = {len(x.fullStructure) for x in self.chordObjects}
		if len(l) != 1:
			print("To use the .sequence method, ChordArray needs each object to have the same number of notes.")
					
		
		for c in self.chordObjects:
			self.array.append([c, copy.deepcopy(c.fullStructure)])
		"""
		
	def sequence(self, transformation = ''):
		
		self.make_array()
		self.mainObject.transform(transformation)
		
		"""TODO: the following doesn't work for direct permutations, that needs to be added"""
		permDict, transpDict = self.mainObject.find_transposition_and_permutation()
		
		for i, data in enumerate(self.array):
			c, theList = data
			if c == self.mainObject: continue
			newList = [None] * len(theList)
			for voiceNumber, noteIndex in enumerate(theList):
				newList[permDict[voiceNumber]] = noteIndex
				if type(noteIndex) is int: 												
					c.currentNotes[noteIndex] += transpDict[voiceNumber]
			"""
			Q: SHOULD ONE UPDATE THE OBJECT'S FULLSTRUCTURE? OR JUST LEAVE IT AS IS?
				leaving it keeps the NHTs intact
				maybe a chord should have a list of all its dependent NHTs, along with whether they update with permutation?
					
			"""
			self.array[i] = [c, newList]
	
	def get_actual_permutation(self, oldList, newList, permDict, transpDict):
		oldList2 = [x for x in oldList if type(x) is not list]
		newList2 = [x for x in newList if type(x) is not list]
		newPermDict = {}
		newTranspDict = {}
		for i, n in enumerate(oldList):
			if type(n) is list: continue
			oldIndex2 = oldList2.index(n)
			transp = transpDict.get(i, 0)
			newIndex2 = newList2.index(n + transp)
			newPermDict[oldIndex2] = newIndex2
			newTranspDict[oldIndex2] = transp
		return newPermDict, newTranspDict
			
	def broadcast(self, transformation = ''):
		"""apply the transformation to every note"""
		for c in self.chordObjects:
			c.transform(transformation)

class VelocityGenerator(Default):	
	defaults = 	{'generators': ['random_choice', 'accent'],
				'additiveGenerators': [],
				'multiplicativeGenerators': [],
				
				'randomRange': [54, 76],
				'accentHisto': {10: 20},
				'accentPct': 10,
				
				'fixedListCounter': -1,
				'fixedList': [],
				'makeList': 0,
								
				'currentFadeFactor': 1.,
				'fadeDelta': .05,
				'fadeOutFunction': None,			# upon completion of fade
				
				'scaleFactor': .8,
				
				'lastWaveVelocity': 36,
				'waveRange': [36, 64],
				'waveDelta': 5,
				'waveDeltaRange': [3, 8],
				'resetWaveBoundaries': True,
				'lowerWaveBoundary': [24, 40],
				'upperWaveBoundary': [76, 96],
				'waveMode': 0						# 0 = reverse, 1 = crec only, -1 = decresc.
				}
	
	def __init__(self, *args, **kwargs):
		"""
		TODO: add a single add_generator function?
		"""
		
		self.set_default_arguments(kwargs)
		
		self.get_generators(self.generators, self.additiveGenerators)
		if self.multiplicativeGenerators:
			self.get_generators(self.multiplicativeGenerators[:], self.multiplicativeGenerators)
		
		"""automatically makes a list of repeating values"""
		if self.makeList:
			self.make_list(self.makeList)
	
	def get_generators(self, gNames, gList):
		gList[:] = []
		
		"""put this in its own routine, so the velocity functions can change on the fly"""
		for funcName in gNames:
			if type(funcName) is str:
				if hasattr(self, funcName):
					gList.append(getattr(self, funcName))
				else:
					print('ARCA ERROR: VelocityGenerator has no attribute', funcName)
			elif type(funcName) is int:
				gList.append(VelocityGenerator.funcgenerator(funcName))
			else:
				gList.append(funcName)
		
	def get_velocity(self):
		self.lastVelocity = 0
		for funcName in self.additiveGenerators:
			self.lastVelocity += funcName()
		for funcName in self.multiplicativeGenerators:
			self.lastVelocity = int(.5 + self.lastVelocity*funcName())
		self.lastVelocity = min(max(0, self.lastVelocity), 127)
		return self.lastVelocity
	
	def random_choice(self):
		return random.randrange(*self.randomRange)
	
	def accent(self):
		if random.randrange(0, 100) < self.accentPct:
			return weighted_choice(self.accentHisto)
		return 0
	
	def make_list(self, listLen = 3):
		self.fixedList = []
		for i in range(listLen):
			self.fixedList.append(self.get_velocity())
		self.get_generators(['fixed_list'], self.additiveGenerators)
	
	def fixed_list(self):
		self.fixedListCounter = (self.fixedListCounter + 1) % len(self.fixedList)
		return self.fixedList[self.fixedListCounter]
		
	def scale(self):
		return self.scaleFactor
	
	def waves(self):
		"""THIS IS KLUDGY; ON FADE OUT, SHOULD ELIMINATE ALL ADDITIVE GENERATORS"""
		if self.currentFadeFactor != 1.:
			return self.lastWaveVelocity
		newVel = self.lastWaveVelocity + self.waveDelta
		self.lastWaveVelocity = newVel
		if newVel > self.waveRange[1]:
			if self.resetWaveBoundaries:
				self.waveRange[0] = random.randrange(*self.lowerWaveBoundary)
			if self.waveMode == 0:
				self.waveDelta = -abs(random.randrange(*self.waveDeltaRange))
			elif self.waveMode == 1:
				self.lastWaveVelocity = self.waveRange[0]
		elif newVel < self.waveRange[0]:
			if self.resetWaveBoundaries:
				self.waveRange[1] = random.randrange(*self.upperWaveBoundary)
			if self.waveMode == 0:
				self.waveDelta = abs(random.randrange(*self.waveDeltaRange))
			elif self.waveMode == -1:
				self.lastWaveVelocity = self.waveRange[1]
		return newVel
	
	def apply_fade(self):
		self.currentFadeFactor = self.currentFadeFactor + self.fadeDelta
		if self.currentFadeFactor <= 0:
			if self.fadeOutFunction:
				self.fadeOutFunction()
			self.currentFadeFactor = 1.
			self.fadeDelta = 0
			self.additiveGenerators = [lambda: 0]
			self.multiplicativeGenerators = []
			return 0
		return self.currentFadeFactor
	
	def constant_function(self):
		return 0
	
	def simple_apply_fade(self):
		return self.lastVelocity - self.simpleFadeDelta
	
	"""the old fade out was creating weird artifacts, this is designed to be really simple"""
	def fade_out_OLD(self, fadeTime = 15, endFunc = None):
		self.simpleFadeDelta = max(1, int(.5 + self.lastVelocity/fadeTime))
		self.additiveGenerators = [self.simple_apply_fade]
		self.multiplicativeGenerators = []
		self.fadeOutFunction = endFunc
	
	"""need a fade_in as well"""
	def fade_out(self, fadeTime = 25, endFunc = None):
		self.fadeDelta = -1./fadeTime
		self.fadeOutFunction = endFunc
		self.currentFadeFactor = 1.
		self.multiplicativeGenerators.append(self.apply_fade)
	
	@staticmethod
	def funcgenerator(f):
		return lambda: f

class RhythmGenerator(Default):
	"""these generators are all multiplicative"""
	defaults = 	{'generators': ['random_choice'],
	
				'totalBeats': -1,
				'augmentation': 0,
				
				'playPct': 50,
				
				'dustFunc': False,		# allows for randomizing the rhythmic pattern
				
				'durationHisto': {1: 50, 2: 25, 3: 5},
				'durationHistoCounter': 0,
				'durationHistoModulus': 2,			# redundant, but here for clarity
				
				'fixedListCounter': -1,
				'fixedList': [],
				'resetList': False,
				'makeList': 0,
								
				'currentFadeFactor': 1.,
				'fadeDelta': .05,
				'fadeOutFunction': None,			# upon completion of fade
				
				'finalRests': 0
				}
	
	def __init__(self, *args, **kwargs):
		"""
		TODO: add a single add_generator function?
		"""
		
		self.set_default_arguments(kwargs)
		
		self.durationHistoModulus = weighted_choice(self.durationHisto)
		self.get_generators(self.generators[:])
		
		"""automatically makes a list of repeating values"""
		if self.makeList:
			self.make_list(self.makeList)
	
	def get_generators(self, gList):
		self.generators = []
		
		"""put this in its own routine, so the velocity functions can change on the fly"""
		for funcName in gList:
			if type(funcName) is str:
				if hasattr(self, funcName):
					self.generators.append(getattr(self, funcName))
				else:
					print('ARCA ERROR: RhythmGenerator has no attribute', funcName)
			else:
				self.generators.append(funcName)
		
	def get_beat(self):
		"""returns 0 if it should play, 1 if not """
		self.totalBeats += 1
		if self.augmentation and self.totalBeats % self.augmentation != 0:
			return 0
		
		totalVel = 1
		for funcName in self.generators:
			totalVel *= funcName()
		if totalVel > 0:
			if self.dustFunc:
				return self.dustFunc()
			else:
				return 1
		return 0
	
	def random_choice(self):
		return int(random.randrange(0, 100) <= self.playPct)
	
	def histogram(self):
		beat = int(self.durationHistoCounter % self.durationHistoModulus == 0)
		self.durationHistoCounter += 1
		if self.durationHistoCounter >= self.durationHistoModulus:
			self.durationHistoModulus = weighted_choice(self.durationHisto)
			self.durationHistoCounter = 0
		return beat
		
	def make_list(self, parameters = 3):
		self.fixedList = []
		if type(parameters) is int:
			for i in range(parameters):
				self.fixedList.append(self.get_beat())
		elif type(parameters) is list:
			self.fixedList = [1]*parameters[0] + [0]*parameters[1]	# parameters is note, rest
			random.shuffle(self.fixedList)
		self.fixedList += [0]*self.finalRests
		self.get_generators(['fixed_list'])
	
	def fixed_list(self):
		self.fixedListCounter = (self.fixedListCounter + 1) % len(self.fixedList)
		if self.fixedListCounter and self.resetList:
			random.shuffle(self.fixedList)
		return self.fixedList[self.fixedListCounter]
	
	def apply_fade(self):
		self.currentFadeFactor = self.currentFadeFactor + self.fadeDelta
		if self.currentFadeFactor <= 0:
			if self.fadeOutFunction:
				self.fadeOutFunction()
			self.currentFadeFactor = 1.
			self.fadeDelta = 0
			return 0
		return int(random(0, 100) <= self.currentFadeFactor)
	
	"""need a fade_in as well"""
	def fade_out(self, fadeTime = 25, endFunc = None):
		self.currentFadeFactor = 100.
		self.fadeDelta = -100./(fadeTime + 1)
		self.fadeOutFunction = endFunc
		self.generators.append(self.apply_fade)



MAXPROGRAMLENGTH = 1000								# to prevent endless loops in Program, should be larger

class Program():
	"""
	Use text to create musical structures
	"""
	defaultObjectArgs = {'Scale': {'useMidi': True}, 'Chord': {'useMidi': True, 'useGroups': True}, 'Arpeggiation': {'useMidi': True, 'useGroups': True}}
	
	def __init__(self, startString = 'currentcode.txt', defaultPath = DEFAULTDIRECTORY, autoRun = True):
		global lastProgram
		
		lastProgram = self
		
		if startString.endswith('.txt'):
			if not startString.count('/'):
				fileName = defaultPath + startString
			with open(fileName, 'r') as myFile:
				startString = myFile.read()
		self.startString = startString
		
		self.objectDict = {}
		
		self.defaultTimeline = Timeline('', self)
		
		self.timeLines = [self.defaultTimeline]
		self.macroDict = {}
		self.masterOutput = {}
		self.noteChanges = {}
		self.groupByTimelines = True
		self.streamOrder = []
		
		self.timeSignature = [4, 4]
		self.beatDuration = .5
		self.measureLength = 8
		self.unitDuration = .5
		self.bps = 2
		
		self.spelling = 0
		
		self.objectsByTimeline = {}
		
		self.tiedNotes = []
		
		self.parse_string(startString)
		
		if not self.streamOrder:
			self.streamOrder = list(range(len(self.timeLines) - 1)[::-1])
		
		self.streamOrder = [0] + [x + 1 for x in self.streamOrder]  					# because of the defaultTimeline	
		
		if autoRun:
			self.output_timelines()
	
	def get_max_time(self, streamObject):
		highestOffset = 0
		maxTime = 0
		data = self.masterOutput[streamObject]
		for offset, chords in data.items():
			for noteData in chords:
				notes, velocities, duration = noteData
				highestOffset = max(offset, highestOffset)
				maxTime = max(maxTime, offset+duration)
		return highestOffset, maxTime
	
	def get_average_pitch(self, streamObject):
		pitchList = []
		data = self.masterOutput[streamObject]
		for offset, chords in data.items():
			for noteData in chords:
				notes, velocities, duration = noteData
				if notes and notes != [None]:
					pitchList += notes
		if not pitchList: return 0
		return sum(pitchList)/len(pitchList)
	
	def show_text(self):
		for obj, data in self.masterOutput.items():
			print('OBJECT', obj)
			for offset, chords in data.items():
				for noteData in chords:
					notes, velocities, duration = noteData
					print('  offset', offset, 'notes', notes, 'velocities', velocities, 'duration', duration)
	
	def show_XML(self, fileName = DEFAULTOUTPUT, autoOpen = True):
		
		"""
		TODO: support key signature
		"""
		
		maxTime = 0
		highestOffset = 0
		
		for obj in self.masterOutput:
			o, t = self.get_max_time(obj)
			maxTime = max(t, maxTime)
		
		for obj in self.masterOutput:
			o, t = self.get_max_time(obj)
			if self.masterOutput[obj] and t < maxTime:
				self.masterOutput[obj][t] = [[None, [64], maxTime - t]]
		
		self.totalParts = []
		writer = DTXML.Writer()
		writer.spelling(direction = self.spelling)
		
		theParts = list(self.masterOutput.items())
		
		for partNumber in self.streamOrder:
			
			kwargs = {'TimeSignature': self.timeSignature}
			
			obj, data = theParts[partNumber]
			newData = {}
		
			if not data:
				self.totalParts.append(False)
				continue
			
			aPitch = self.get_average_pitch(obj)
			
			if aPitch < 60:
				kwargs['Clef'] = ['F', 4]
			
			for offset, chords in data.items():
				for noteData in chords:
					notes, velocities, duration = noteData
					velocity = velocities[0]
					newData[offset*self.beatDuration] = notes, duration*self.beatDuration
			writer.add_part(newData, None, **kwargs)
		
		if fileName is False:	
			writer.write()
		elif fileName is None:
			return writer.get_XML_as_string()
		else:
			writer.write(fileName)
			
		if writer.fileName and autoOpen:
			import subprocess
			subprocess.call(('open', writer.fileName))
	
	def show(self, outputFormat = '', **kwargs):
		outputFormat = outputFormat.upper().replace(' ', '')
		if not outputFormat or outputFormat.startswith('X'):
			self.show_XML(**kwargs)
			return
		else:
			if outputFormat.startswith('T'):
				self.show_text()
				return
			elif outputFormat.startswith('M'):
				self.show_music21()
				return
			
	def show_music21(self):
		try:
			import music21							# TODO: enclose in a try loop, provide another way of showing?
		except:
			print("MUSIC21 IMPORT ERROR")
			return False
		
		self.totalParts = []
		
		for obj, data in self.masterOutput.items():
		
			if not data:
				self.totalParts.append(False)
				continue
		
			p = music21.stream.Part()
			p.insert(0, music21.meter.TimeSignature(f'{self.timeSignature[0]}/{self.timeSignature[1]}'))
			self.totalParts.append(p)
		
			for offset, chords in data.items():
				for noteData in chords:
					notes, velocities, duration = noteData
					velocity = velocities[0]
					if notes is None or notes == [None]:						# todo: how does brython deal with this
						c = music21.note.Rest(quarterLength = duration*self.beatDuration)
					else:
						c = music21.chord.Chord(notes, quarterLength = duration*self.beatDuration)
					p.insert(offset*self.beatDuration, c)
				
		self.output = music21.stream.Stream()
		maxLength = max([p.quarterLength for p in self.totalParts if p])
	
		for i in self.streamOrder:
			p = self.totalParts[i]
			if p:
				self.output.insert(0, p)
				d = maxLength - p.quarterLength
				if d > 0:
					p.append(music21.note.Rest(quarterlength = d))
				p.makeVoices(inPlace = True)
		
		self.output.show()
	
	def consolidate_output(self):
		self.consolidatedOutput = {}
		for obj, data in self.masterOutput.items():
			if not data: continue
			for offset, chords in data.items():
				for noteData in chords:
					notes, velocities, duration = noteData
					l = len(velocities)
					if notes is None:
						continue
					for i, n in enumerate(notes):
						if n is None: continue
						nData = [n, velocities[i % l], duration]
						self.consolidatedOutput.setdefault(offset, []).append(nData)
		return self.consolidatedOutput
						
	def objects(self):
		import __main__
		for k, val in self.objectDict.items():
			setattr(__main__, k, val)
	
	def regularize_offsets(self):
		newDict = {}
		for obj, data in self.masterOutput.items():
			newDict[obj] = {}
			for offset, chords in data.items():
				if type(offset) is float:
					o = round(offset, 2)
					offset = int(o) if o.is_integer() else o
				for c in chords:
					newDict[obj].setdefault(offset, []).append(c)
		self.masterOutput = newDict
		self.resolve_ties()
	
	def resolve_ties(self):
		for offset, obj, note, dur in self.tiedNotes:
			originalChordList = self.masterOutput[obj][offset]
			originalNote = False
			for noteData in originalChordList:
				if note in noteData[0]:
					originalNote = noteData
					break
			if not originalNote:
				print("ARCA ERROR: can't find original tied note!")
				continue
			target = self.masterOutput.get(obj, {}).get(offset + dur, [])
			targetNote = False
			for noteData in target:			# notes, velocities, duration = noteData
				if note in noteData[0]:
					targetNote = noteData
					break
			if not targetNote:
				print("ARCA ERROR: can't find destination tied note!")
				continue
			originalNote[2] += targetNote[2]
			targetNote[0].remove(note)
			if not targetNote[0]:
				target.remove(targetNote)
	
	def play_minimal_offset(self):
		global minOffset
		
		minElements = []
		for x in self.timeLines:
			if x.actionItems:
				minElements.append([min(x.actionItems), x])
		
		if not minElements:
			return False
		
		minOffset = min(minElements, key = lambda x: x[0])[0]
		
		timeLinesToPlay = [x[1] for x in minElements if x[0] == minOffset]
		
		if not timeLinesToPlay:
			print("ARCA ERROR: I SHOULDN'T BE HERE (play_minimal_offset)!!!")
			return False
		
		for t in timeLinesToPlay:
			t.play_events_at_offset(minOffset)
			self.objectsByTimeline[t.currentObject] = self.timeLines.index(t)
			
		return True
	
	def output_timelines(self):
		
		going = True
		counter = 0
		self.masterOutput = {o:{} for o in self.timeLines}
		
		while going and counter < MAXPROGRAMLENGTH:
			going = self.play_minimal_offset()
			counter += 1
		
		self.regularize_offsets()
		self.add_and_subtract_notes()
		"""TODO: add rests to make all notes the same length"""
	
	def add_and_subtract_notes(self):
		
		for o, changesToMake in self.noteChanges.items():
			addedNotes, subtractedNotes = changesToMake
			for obj, noteData in addedNotes:
				chordList = self.masterOutput[obj].setdefault(o, [])
				if not chordList:
					chordList.append(noteData)
				else:
					chordList[0][0] += noteData[0]
			for obj, noteData in subtractedNotes:
				chordList = self.masterOutput[obj].setdefault(o, [])
				"""TODO: iterate over all chords at this timepoint"""
				for n in noteData:
					if n in chordList[0][0]:
						chordList[0][0].remove(n)
	
	def parse_string(self, inputString = ''):
		
		"""Split into lines.  
			Continuation lines begin with a tab"""
		inputString = inputString.replace('\n\t', ' ').replace('\n', ';').replace('\t', ' ').replace('&', '!next ').replace('!', '! ').replace('\r', '')
		inputString = self.purify_string(inputString)
		inputSplit = inputString.split(';')
		
		"""Filter comments"""
		for i, s in enumerate(inputSplit):
			if s.count('//'):
				inputSplit[i] = s[:s.index('//')]
		
		for s in inputSplit:
			"""TODO: allow for time signature = beat duration = etc. (two words)"""
			sSplit = s.split()
			if not sSplit: continue
			if sSplit[0].startswith('@'):
				if len(sSplit) == 1: continue
				targetBeat = parse_number(sSplit[0][1:])
				sSplit = sSplit[1:]
			else:
				targetBeat = 0
			if len(sSplit) > 1 and sSplit[1] == '=':
				result = self.parse_assignment(sSplit)
				if result:
					self.defaultTimeline.add_event(targetBeat, result)
			elif sSplit[0].endswith('Stream:'):
				self.timeLines.append(Timeline(sSplit, self))
	
	def update_tempo(self):
		"""need to call this every time the tempo changes because we want up-to-date parsing of timelines"""
		self.beatLength = 4/self.timeSignature[1]											# quarters/beatSubdivision
		unitsPerQuarter = 1/self.beatDuration												# programUnits/quarter (TODO: rename beatDuration "divisions"?)
		self.measureLength = unitsPerQuarter*self.beatLength*self.timeSignature[0]			# masure length in program units
		if self.timeSignature[0] != 3 and self.timeSignature[0] % 3 == 0:					# 3/8, 3/4, etc. ar 								
			self.beatLength = 3 * self.beatLength											# quarters/beat
		elif self.timeSignature[1] % 5 == 0:
			self.beatLength = 3 * self.beatLength											# quarters/beat
		self.unitDuration = 1/(unitsPerQuarter*self.beatLength*self.bps)					# 1/(units/quarter * quarters/beat * beats/second)
	
	def parse_assignment(self, sSplit):
		sSplit = [self.apply_macro(x) for x in sSplit]
		varName = sSplit[0]
		remainder = ' '.join(sSplit[2:])
		varUp = varName.upper()
		if varUp == 'TEMPO':
			self.bps = parse_number(sSplit[2])/60											# beats/second
			self.update_tempo()
			return
		elif varUp == 'TIMESIGNATURE':
			numerator, denominator, *junk = sSplit[2].split('/')
			self.timeSignature = [int(numerator), int(denominator)]
			self.update_tempo()
			#print(self.timeSignature, f'Measure length is {self.measureLength}')
			return
		elif varUp == 'BASICDURATION' or varUp == 'BEATDURATION':
			self.beatDuration = parse_number(sSplit[2])
			self.update_tempo()
			return
		elif varUp == 'STREAMORDER':
			self.streamOrder = [parse_number(x) for x in sSplit[2:]]
			return
		elif varUp == 'SPELLING':
			direction  = sSplit[2].upper()
			if direction.count('SHARP'):
				self.spelling = 1
			elif direction.count('FLAT'):
				self.spelling = -1
			return
		#print(sSplit, remainder)
		firstParen = remainder.index('(')
		lastParen = [i for i in range(len(remainder)) if remainder[i] == ')'][-1]
		objectName = remainder[:firstParen]
		
		if objectName == 'Macro':
			
			self.macroDict[varName] = remainder[firstParen+1:lastParen]
			return False
		
		args = [varName, objectName]
		
		l, d = self.parse_arguments(remainder[firstParen+1:lastParen])
		
		for abbrev, fullKeyword in [['chain', 'chainMotive']]:
			if abbrev in l:
				l.remove(abbrev)
				d[fullKeyword] = True
		
		args = [x for x in args + l if x]
		kwargs = copy.deepcopy(Program.defaultObjectArgs.get(objectName, {}))
		for k, val in d.items():
			kwargs[k] = val

		return([self.defaultTimeline.make_object, args, kwargs])
	
	def parse_arguments(self, s):
		args = []
		kwargs = {}
		parameters = parse_parentheses(s.split(','), ',')
		
		for p in parameters:
			"""self.apply_macro here???"""
			p = p.rstrip(' ').lstrip(' ')
			if p.count('='):
				k, val, *junk = p.split('=')
				k = k.replace(' ', '')
				if k in ['container', 'chord']:
					val = val.replace(' ', '')
					val = self.objectDict.get(val, val)
				elif k in ['durations', 'length']:
					val = [parse_number(x) for x in val.split()]
				elif k in ['progression']:
					val = parse_parentheses(val.split())
				kwargs[k] = val
			else:
				args.append(p)
		return args, kwargs
		
	def get_object(self, objectOrName):
		if type(objectOrName) is str:
			objectOrName = self.objectDict.get(objectOrName)
			if not objectOrName:
				print("ARCA ERROR: Can't find object with name", objectOrName)
		return objectOrName
			
	def apply_macro(self, s):
		for k, val in self.macroDict.items():
			s = s.replace(k, val)
		return s
	
	def purify_string(self, s):
		return ' '.join(s.split()).replace(' )', ')').replace('( ', '(').replace('\r', '').rstrip().lstrip()

class Timeline(Default):
	
	defaults = 	{	'currentOffset': 0,
					'currentObject': False,
					'program': False,
					'actionItems': {},
					'targetBeat': 0
				}
				
	"""
	
	Timelines are used by the Program class.  They are essentially a dictionary whose keys are times and whose values are lists of function calls.  
	What is complex is that you don't necessarily know at what time events should happen: with a command like !!! you know an object should output three
	times but not how long the object is, since it can change at runtime.
				
	The main structure is self.actionItems, which has the format 
		offset: listOfSequences
	
	A sequence is itself a list of commands to be executed sequentially, like !!!! 
	
	A command has the format:
		[function, args, kwargs]
	
	So the minimal timeLine is:		
				offset: [[[function, args, kwargs]]]
				
	"""
				
	def __init__(self, *args, **kwargs):
		
		self.set_default_arguments(kwargs)
		
		if args and len(args) >= 2:
			self.startCommands = args[0]
			self.program = args[1]
			if self.startCommands:
				self.parse_command_list(self.startCommands, self.targetBeat)
			
	def add_sequence(self, offset, sequence):						# a sequence is a list of function calls
		self.actionItems.setdefault(offset, []).append(sequence)
	
	def add_event(self, offset, event):								# an event is a single function call, so it needs to be embedded in a list
		self.actionItems.setdefault(offset, []).append([event])
		
	def play_events_at_offset(self, offset):
		
		if offset not in self.actionItems:
			return False
		
		sequenceList = self.actionItems.get(offset, [])
		
		"""TODO: fix Brython kludge"""
		try:
			del self.actionItems[offset]
		except:
			i = int (offset)
			del self.actionItems[i]			# if i in self.actionItems ... [???]
		
		for sequenceNumber, sequence in enumerate(sequenceList):
			self.currentOffset = offset
			for i, action in enumerate(sequence):
				func, args, *kwargs = action
				kwargs = kwargs[0] if kwargs else {}
				func(*args, **kwargs)
				if self.currentOffset != offset:
					break
			if self.currentOffset != offset and len(sequence) > i+1:
				self.add_sequence(self.currentOffset, sequence[i+1:])
		
		"""try:
			
		except:
			pass"""
		
	def parse_command_list(self, inputList, targetBeat = 0):
		
		"""Assemble a list of function calls from text"""
		
		sSplit = parse_parentheses(inputList)
		
		out = []
		
		for item in sSplit:
			item = item.rstrip().lstrip()
			if item.startswith('('):
				out[-1] = out[-1] + item
			else:
				out.append(item)
		
		sSplit[0] = sSplit[0].replace('Stream:', ':')
		if all([x.isdigit() for x in sSplit[0]]):
			"""TODO: add code here to name streams in a way that orders them"""
			pass
			
		sequence = []
		self.add_sequence(targetBeat, sequence)
		
		for rawItem in sSplit:
			
			itemList = self.program.apply_macro(rawItem)
			
			if itemList.startswith('{'):
				x = self.program.parse_assignment(itemList[1:-1].split())
				if x:
					sequence.append(x)
				continue
			
			if itemList == rawItem:
				itemList = [itemList]
			else:
				itemList = parse_parentheses(itemList.split())
			
			for item in itemList:
			
				if item.endswith(':'):
					sequence.append([self.change_object, [item[:-1]]])
						
				elif item.startswith('@+') or item.startswith('@-'):
					sequence.append([self.change_beat, [parse_number(item[1:])]])
				
				elif item.startswith('@m+') or item.startswith('@m-'):
					targetBeat = (parse_number(item[2:]) - 1)*self.program.measureLength
					if type(targetBeat) is float and targetBeat.is_integer():
						targetBeat = int(targetBeat)
					sequence.append([self.change_beat, [targetBeat]])
				
				elif item.startswith('@m'):
					item = item[2:]
					if item.count('+'):
						i = item.index('+')
						addedBeats = parse_number(item[i:])
						item = item[:i]
					elif item.count('-'):
						i = item.index('-')
						addedBeats = parse_number(item[i:])
						item= item[:i]
					else:
						addedBeats = 0
					targetBeat = (parse_number(item) - 1)*self.program.measureLength + addedBeats
					if type(targetBeat) is float and targetBeat.is_integer():
						targetBeat = int(targetBeat)
					sequence = []
					self.add_sequence(targetBeat, sequence)
				
				elif item.startswith('@'):
					targetBeat = parse_number(item[1:])
					sequence = []
					self.add_sequence(targetBeat, sequence)
					
				elif item.endswith('!'):
					if item.count('['):
						indices = ''.join([x for x in item if x.isdigit() or x in ':-']).split(':')
						indices = [int(x) if x else None for x in indices]
						if len(indices) == 1:
							indices = [indices[0], indices[0] + 1, 1]
						elif len(indices) == 2:
							indices.append(None)
						if indices[0] is None:
							indices[0] = 0
						"""can't correct indices[1] until we have the notes"""
						if indices[2] is None:
							indices[2] = 1
					else:
						indices = None
					sequence.append([self.play_object, [indices]])
						
				elif item.startswith('+'):
					if item[1] == '(':
						iSplit = item.split(',')
						notes = [parse_lettername(x) for x in iSplit[0][2:].split()]
						pList = ','.join(iSplit[1:]) if len(iSplit) > 1 else []
					else:
						iSplit = item.split('(')
						notes = [parse_lettername(iSplit[0][1:])]
						pList = iSplit[1] if len(iSplit) > 1 else []
					d = {}
					if pList:
						l, d = self.program.parse_arguments(pList.replace(')', ''))
						for abbrev, fullKeyword in [['tied', 'tied']]:
							if abbrev in l:
								l.remove(abbrev)
								d[fullKeyword] = True
					length = d.get('length', 0)
					if type(length) is list:
						length = length[0]
						d['extraDurations'] = length
					if d.get('length'):
						del d['length']
					if length:
						sequence.append([self.add_notes, [notes, length], d])
					else:
						"""TODO: a different function that takes the chord's duration"""
						sequence.append([self.add_notes, [notes, None], d])
				
				elif item.startswith('-'):
					if item[1] == '(':
						iSplit = item.split(',')
						notes = [parse_lettername(x) for x in iSplit[0][2:].split()]
						pList = ','.join(iSplit[1:]) if len(iSplit) > 1 else []
					else:
						iSplit = item.split('(')
						notes = [parse_lettername(iSplit[0][1:])]
						pList = iSplit[1] if len(iSplit) > 1 else []
					sequence.append([self.remove_notes, [notes]])
			
				elif item.startswith('newscale'):
					params = item.split('(')[1].rstrip(') ')
					params = [parse_lettername(x) for x in params.split()]
					if type(params) is list:
						if len(params) == 1:
							params = params[0]
						else:
							params = sorted(params)
					sequence.append([self.newscale, [params]])
				
				elif item.startswith('addnote'):
					# ADD NOTE TO CURRENT CHORD
					params = item.split('(')[1].rstrip(') ')
					notes = [parse_lettername(x) for x in params.split()]
					sequence.append([self.add_note_to_chord, [notes]])
				
				elif item.startswith('removenote'):
					params = item.split('(')[1].rstrip(') ')
					notes = [parse_lettername(x) for x in params.split()]
					sequence.append([self.remove_note_from_chord, [notes]])
				
				elif item.startswith('next'):
					sequence.append([self.apply_progression, []])
				
				elif item.startswith('changeduration'):
					params = []
					for pair in item.split(','):
						item1, item2, *junk = [parse_number(x) for x in item.split()]
						params.append([int(item1), item2])
					
					sequence.append([self.change_duration, [params]])
				
				else:
					sequence.append([self.transform, [item]])
	
	def apply_progression(self):
		self.currentObject.next()
	
	def change_object(self, objectName):
		newObject = self.program.objectDict.get(objectName)
		if not newObject:
			print("ARCA ERROR: Can't find object named", objectName)
			return False
		self.currentObject = newObject
	
	def change_beat(self, timeDelta = 0):
		self.currentOffset += timeDelta
	
	def play_object(self, indices = None):
		
		"""
		This is kind of a sloppy hack
		
		TODO:
			1. I have not really tested the situation where you want to play only, e.g. every other note in a motive (might work as is)
			2. It would be better to build the durations note by note, to allow the motive durations to change during the motive's lifespan; 
				right now the durations are fixed at launch.  This is a relatively simple matter of just adding the playnote events one by one.
		
		"""
		durations = self.currentObject.get_durations()
		
		if indices:
			if indices[1] is None:
				indices[1] = len(durations)
			#if indices[0]:						# ELIMINATED THIS.  MIGHT CAUSE PROBLEMS???
			functionToCall = self.playnote_incomplete
			iterationObj = range(indices[0], indices[1], indices[2])
		if indices is None:
			iterationObj = range(len(durations))
			functionToCall = self.playnote
		
		for i in iterationObj: #, dur in enumerate(durations):
			dur = durations[i]
			self.add_event(self.currentOffset, [functionToCall, [i, dur]])
			self.currentOffset += dur
			functionToCall = self.playnote
	
	def playnote(self, index, dur):
		obj = self.currentObject
		theNotes = obj.get_group(index)
		velocities = [64]							# placeholder
		noteData = [theNotes, velocities, dur]
		self.program.masterOutput.setdefault(self, {}).setdefault(self.currentOffset, []).append(noteData)
	
	"""TODO: I don't think we need this if get_group can backfill, calculating all notes up to the desired one"""
	def playnote_incomplete(self, index, dur):
		obj = self.currentObject
		for i in range(0, index + 1):
			theNotes = obj.get_group(i)
		velocities = [64]							# placeholder
		noteData = [theNotes, velocities, dur]
		self.program.masterOutput.setdefault(self, {}).setdefault(self.currentOffset, []).append(noteData)
		
	def transform(self, cmd):
		self.currentObject.transform(cmd)
	
	def newscale(self, theScale, **kwargs):
		
		scaleObject = self.currentObject
		if not scaleObject: return False
		
		for name, o in self.program.objectDict.items():
			if o.container == scaleObject and type(o) != Arpeggiation:
				o.get_MIDI()
		
		scaleObject.change_scale(theScale, container = scaleObject.container, **Program.defaultObjectArgs['Scale'])  	# more parameters needed?
		newModulus = len(scaleObject)
		
		for name, o in self.program.objectDict.items():
			if o.container == scaleObject and type(o) != Arpeggiation:
				o.restore_from_MIDI()
				o.get_MIDI()
				o.modulus = newModulus
				if type(o) is Chord:
					for a in o.dependentObjects:
						a.modulus = newModulus
						
	def change_duration(self, params = []):
		obj = self.currentObject
		l = len(obj.durations)
		for durationNumber, newDuration in params:
			obj.durations[durationNumber] = newDuration
			
	def add_note_to_chord(self, notes):
		
		obj = self.currentObject
		if type(obj) == Chord:
			if obj.container:
				notes = [obj.container.scale_degree(x) for x in notes]
			obj.currentNotes += notes
			obj.update()
		"""TODO: need code for scale here"""
	
	def add_notes(self, theNotes, length = None, **kwargs):
		
		if length is None:
			dur = 1
		else:
			dur = length
			
		obj = self.currentObject
		
		if type(theNotes) is str:
			theNotes = parse_lettername(theNotes)
		if type(theNotes) is not list:
			theNotes = [theNotes] 
			
		velocities = [64]							# placeholder
		noteData = [theNotes, velocities, dur]
		if kwargs.get('tied'):
			for note in theNotes:
				self.program.tiedNotes.append([self.currentOffset, self, note, dur])
		
		if length: 
			self.program.masterOutput.setdefault(self, {}).setdefault(self.currentOffset, []).append(noteData)
			self.currentOffset += length
		else:
			self.program.noteChanges.setdefault(self.currentOffset, [[], []])[0].append([self, noteData])
	
	def remove_notes(self, theNotes, **kwargs):
		self.program.noteChanges.setdefault(self.currentOffset, [[], []])[1].append([self, theNotes])
		pass
	
	def make_object(self, varName, objectName, *args, **kwargs):
		global objectDict
		if objectName in ['Scale', 'Chord', 'Arpeggiation']:
			container = kwargs.get('container')
			if type(container) is str:
				kwargs['container'] = self.program.get_object(container)
			
			chord = kwargs.get('chord')
			if type(chord) is str:
				kwargs['chord'] = self.program.get_object(chord)
			
			kwargs['objectDict'] = self.program.objectDict
			self.program.objectDict[varName] = globals()[objectName](*args, **kwargs)
			objectDict = self.program.objectDict
		else:
			print("ARCA ERROR: Can't find object with name", varName)
