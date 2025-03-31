import copy

# a very simple musicXML writer 

DEFAULTOUTPUT = '/Users/denis/Desktop/output1.musicxml'

header = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n<!DOCTYPE score-partwise PUBLIC "-//Recordare//DTD MusicXML 3.0 Partwise//EN" "http://www.musicxml.org/dtds/partwise.dtd">\n<score-partwise version="3.0">'

testDict = {0: ['G4', 1], 1: ['A4', 2], 3: ['B4', 2]}
testDict2 = {0: [67, 1], 1: [69, 2], 3: [71, 2]}

class Writer():
	
	noteDurations = [[0.015625, '256th'], [0.03125, '128th'], [0.0625, '64th'], [0.125, '32nd'], [0.25, '16th'], [0.5, 'eighth'], [1, 'quarter'], [2, 'half'], [4, 'whole'], [8, 'breve'], [16, 'long']]
	accidentals = {-3: 'triple-flat', -2: 'flat-flat', -1: 'flat', 1: 'sharp', 2: 'double-sharp', 3: 'triple-sharp'} #add quarter-tone accidentals
	sharpDict = {0: 'C', 1: 'C#', 2: 'D', 3: 'D#', 4: 'E', 5: 'F', 6: 'F#', 7: 'G', 8: 'G#', 9: 'A', 10: 'A#', 11: 'B'}
	flatDict = {0: 'C', 1: 'Db', 2: 'D', 3: 'Eb', 4: 'E', 5: 'F', 6: 'Gb', 7: 'G', 8: 'Ab', 9: 'A', 10: 'Bb', 11: 'B'}
	
	def __init__(self, fileName = False):
		self.indentLevel = 0
		self.indentString = '  '
		self.divisions = 1440					# divisions of the quarter note
		self.fileName = fileName
		self.lastMeasureNumber = 0
		self.lastPartNumber = 0
		self.lastTimeSignature = [4, 4]
		self.get_measure_length()
		self.parts = {}
		self.partArguments = {}
		self.globalMeasureAttributes = {}
		self.writeString = ''
		self.spellingDict = Writer.sharpDict
	
	def spelling(self, spellingDict = {}, direction = 0):
		if direction > 0:
			self.spellingDict = Writer.sharpDict
		elif direction < 0:
			self.spellingDict = Writer.flatDict
		elif spellingDict:
			self.spellingDict = spellingDict
	
	def make_part_list(self):
		self.output('<part-list>', 1)
		for i, partName in enumerate(self.parts):
			self.output(f'<score-part id="P{i}">', 1)
			self.output(f'<part-name>{partName}</part-name>', 1)
			self.output('</score-part>', -1, -1)
		self.output('</part-list>')
		self.comment_divisor()
	
	def add_part(self, partData, name = None, **kwargs):
		if name is None:
			name = str(len(self.parts))
		self.parts[name] = partData
		self.partArguments[name] = kwargs
		
	def get_XML_as_string(self):
		self.indentLevel = 0
		self.writeString = ''
		self.fileName = ''
		self.file = False
		self.output(header)
		self.make_part_list()
		for name, p in self.parts.items():
			self.write_part(p)
		self.comment_divisor()
		self.output('</score-partwise>', -1)
		return self.writeString
	
	def write(self, fileName = DEFAULTOUTPUT):
		self.indentLevel = 0
		self.fileName = fileName
		if fileName:
			self.file = open(fileName, 'w+')
		else:
			self.file = False
		self.output(header)
		self.make_part_list()
		for name, p in self.parts.items():
			self.write_part(p)
		self.comment_divisor()
		self.output('</score-partwise>', -1)
		if fileName:
			self.file.close()
	
	def write_part(self, partData):
		if not partData: return
		idNumber = list(self.parts.values()).index(partData)
		name = [x[0] for x in self.parts.items() if x[1] == partData][0]
		self.lastMeasureNumber = 0
		kwargs = self.partArguments[name]
		self.output(f'<part id="P{idNumber}">', 0, 1)
		working = partData
		while working:
			working = self.make_measure(working, **kwargs)
			kwargs = {}
		self.output(f'</part>', -1)
	
	def comment_divisor(self):
		self.output('<!--=========================================================-->')
	
	def indent(self, s, changeIndentBefore = 0, changeIndentAfter = 0):
		self.indentLevel += changeIndentBefore
		s = self.indentString * self.indentLevel + s
		self.indentLevel += changeIndentAfter
		return s
	
	def output(self, *args, **kwargs):
		t = self.indent(*args, **kwargs) + '\n'
		if self.file:
			self.file.write(t)
		else:
			self.writeString += t
	
	def make_measure(self, offsetDict, **kwargs):
		self.lastMeasureNumber += 1
		theArgs = {}
		for targetDict in [self.globalMeasureAttributes.get(self.lastMeasureNumber, {}), kwargs]:
			for key, val in targetDict.items():
				theArgs[key] = val
		self.output(f'<measure number="{self.lastMeasureNumber}">')
		self.output('<attributes>', 1)
		self.output(f'<divisions>{self.divisions}</divisions>', 1)
		obj = theArgs.get('KeySignature', None)
		if obj is not None:
			self.output('<key>')
			self.output(f'<fifths>{obj[0]}</fifths>', 1)
			self.output(f'<mode>{obj[1]}</mode>')
			self.output('</key>', -1)
		obj = theArgs.get('TimeSignature', None)
		if (not obj) and self.lastMeasureNumber == 0:
			obj = self.lastTimeSignature
		if obj is not None:
			self.output('<time>')
			self.output(f'<beats>{obj[0]}</beats>', 1)
			self.output(f'<beat-type>{obj[1]}</beat-type>')
			self.output('</time>', -1)
			self.lastTimeSignature = obj
			self.get_measure_length()
		obj = theArgs.get('Clef', None)
		if obj is not None:
			self.output('<clef>')
			self.output(f'<sign>{obj[0]}</sign>', 1)
			self.output(f'<line>{obj[1]}</line>')
			self.output('</clef>', -1)
		self.output('</attributes>', -1)
		
		measureOverflow = False
		totalDuration = 0
		
		for o, data in offsetDict.items():
			if o >= self.measureLength:
				measureOverflow = [True, []]
				break
			if o > totalDuration:
				self.single_pitch(False, o - totalDuration)
			totalDuration = o
			lastOffset = o
			notes, duration, *extra = data											# extra is a dict containing keyword items
			extra = extra[0] if extra else {}
			if duration + totalDuration > self.measureLength:
				d1 = self.measureLength - totalDuration
				d2 = duration + totalDuration - self.measureLength
				measureOverflow = [True, [notes, d2, copy.deepcopy(extra)]]
				extra['tie'] = 'start'
				duration = d1
			self.make_note(notes, duration, **extra)
			totalDuration += duration
			if measureOverflow: break
		
		
		if totalDuration < self.measureLength:
			self.single_pitch(False, self.measureLength - totalDuration)
			
		self.output('</measure>', -1)
		
		newOffsetDict = {}
		if measureOverflow:
			if measureOverflow[1]:
				measureOverflow[1][-1]['tie'] = 'stop'
				newOffsetDict[0] = measureOverflow[1]
			for k, val in offsetDict.items():
				if k >= self.measureLength:
					newOffsetDict[k - self.measureLength] = val
		
		return newOffsetDict
		
	def make_note(self, n = False, duration = 1, **kwargs):
		if n and type(n) is list:
			for i, individualNote in enumerate(n):
				self.single_pitch(individualNote, duration, chord = bool(i != 0), **kwargs)
		else:
			self.single_pitch(n, duration, **kwargs)
		
	def single_pitch(self, n = False, duration = 1, **kwargs):
		"""
		chord = False, tied = True
		"""
		self.output('<note>')
		alter = False
		if n:
			if type(n) is int:
				n = self.spell(n)
			step, alter, octave = parse_letter(n)
			if kwargs.get('chord'): self.output('<chord/>')	
			self.output('<pitch>', 1)	
			self.output(f'<step>{step}</step>', 1)
			if alter:
				self.output(f'<alter>{alter}</alter>')
			self.output(f'<octave>{octave}</octave>')
			self.output('</pitch>', -1)
		else:
			self.output('<rest/>', 1)
		durationInDivisions = duration * self.divisions
		if type(durationInDivisions) is float:
			if durationInDivisions.is_integer():
				durationInDivisions = int(durationInDivisions)
			else:
				print("NON INTEGRAL DURATION IN DIVISIONS HERE, part", self.lastPartNumber, "measure", self.lastMeasureNumber)
				durationInDivisions = int(.5 + durationInDivisions)
		self.output(f'<duration>{durationInDivisions}</duration>')
		tie = kwargs.get('tie')
		if tie:
			self.output(f'<tie type="{tie}"/>')
		self.output(f'<voice>1</voice>')						# not essential
		self.convert_duration_to_XML(duration)
		if n and alter:
			self.convert_accidental_to_XML(alter)
		if tie:
			self.output('<notations>')
			self.output(f'<tied type="{tie}"/>', 1)
			self.output('</notations>', -1)
		self.output(f'</note>', -1)
		
	def convert_accidental_to_XML(self, alter = 0):
		if not alter: return
		s = Writer.accidentals.get(alter)
		if not s:
			print("XML WRITER ERROR, unknown accidental", alter)
		self.output(f'<accidental>{s}</accidental>')
	
	def convert_duration_to_XML(self, duration):
		numDots = 0
		roundLevel = 4
		for i, data in enumerate(Writer.noteDurations):
			if data[0] > duration:
				break
		dur, symbol = Writer.noteDurations[i-1]
		if round(duration, roundLevel) == round(dur, roundLevel):
			numDots = 0
		elif round(duration, roundLevel) == round(dur*1.5, roundLevel):
			numDots = 1
		elif round(duration, roundLevel) == round(dur*1.75, roundLevel):
			numDots = 2
		elif round(duration, roundLevel) == round(dur*1.875, roundLevel):
			numDots = 3
		self.output(f'<type>{symbol}</type>')
		for i in range(numDots):
			self.output('<dot/>')
		
	def get_measure_length(self):
		self.beatLength = 4/self.lastTimeSignature[1]											# quarters/beat unit
		self.measureLength = self.beatLength*self.lastTimeSignature[0]							# masure length in quarter notes
		
	def spell(self, i):
		return f'{self.spellingDict.get(i%12)}{int(i/12) - 1}'
		
def parse_letter(s):
	letter = s[0]
	remainder = s[1:]
	octave = int(''.join([x for x in remainder if x.isdigit()]))
	remainder = remainder.upper()
	alter = remainder.count('#') - remainder.count('B') - remainder.count('-')
	return letter, alter, octave