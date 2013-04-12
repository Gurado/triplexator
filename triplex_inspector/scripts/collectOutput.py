#!/usr/bin/python

import sys, os, traceback
from optparse import OptionParser
import fileinput
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import multiprocessing, Queue

global options, args
global records, loiBed, offBed, targets, tfotargets, annotations

class Job:
	def __init__(self, targetid, tfo_region, tfo_error, lflank, rflank, tsize, tstrand):
		self.targetid = targetid
		self.tfo_region = tfo_region
		self.tfo_error = tfo_error
		self.lflank = lflank
		self.rflank = rflank
		self.tsize = tsize
		self.tstrand = tstrand

class Worker(multiprocessing.Process):
	
	def __init__(self, work_queue, result_queue, chromatin, chromatinFormat):
	
		# base class initialization
		multiprocessing.Process.__init__(self)
		
		# job management stuff
		self.work_queue = work_queue
		self.result_queue = result_queue
		self.kill_received = False
		
		# chromatin data
		self.chromatin = chromatin
		self.chromatinFormat = chromatinFormat
		self.samfile = ""
		self.bigwig = ""
		
		# check if chromatin file has been supplied and create proper filehandle
		if (self.chromatin != "" and self.chromatinFormat == "bam"):
			import pysam # lazy import of required module
			# try open the bam file 
			try:
				self.samfile = pysam.Samfile(self.chromatin, "rb" )
			except:
				exit(1)
		elif (self.chromatin != "" and self.chromatinFormat == "bigwig"):
			from bx.intervals.io import GenomicIntervalReader # lazy import of required module
			from bx.bbi.bigwig_file import BigWigFile  # lazy import of required module
			# try open the bigwig file
			try:
				self.bigwig = BigWigFile(open(self.chromatin))
			except:
				exit(1)

	def __del__(self):
		if (self.chromatin != "" and self.chromatinFormat == "bam"):
			try:
				self.samfile.close()
			except:
				exit(1)	
	
	def run(self):
		while not self.kill_received:
			try:
				# get a task
				job = self.work_queue.get()
				off_target = {}
				processOfftarget(off_target, job.targetid, job.tfo_region, job.tfo_error, job.lflank, job.rflank, job.tsize, job.tstrand, self.samfile, self.bigwig)
				
				# store the result
				self.result_queue.put(off_target)
				self.work_queue.task_done()
			
			except Queue.Empty:
				print 'TIMEOUT' 

def readdata():
	global records, loiBed, offBed, targets, tfotargets, annotations
	
	# read fasta if supplied
	records = {}
	if (options.fastaFile != ""):
		handle = open(options.fastaFile, "rU")
		records = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
		handle.close()	
		if (options.verbose):
			print >> sys.stderr, "FASTA: read %d entries"  % (len(records))
			
	# read region of interest BED file
	loiBed = {}
	try:
		for line in fileinput.input([args[0]]):
			if (len(line.strip())==0 or line[0]=="#"):
				continue
			cols = line.strip().split("\t")
			loiBed[cols[3]] = cols
		if (options.verbose):
			print >> sys.stderr, "LOI BED: read %d entries"  % (len(loiBed))
	except:
		print >> sys.stderr, "[ERROR] LOI BED could not be read : %s" % (args[0])
		exit(1)
		
	# read alternative regions to consider BED file
	offBed = {}
	try:
		for line in fileinput.input([args[1]]):
			if (len(line.strip())==0 or line[0]=="#"):
				continue
			cols = line.strip().split("\t")
			offBed[cols[3]] = cols
		if (options.verbose):
			print >> sys.stderr, "Off BED: read %d entries"  % (len(offBed))
	except:
		print >> sys.stderr, "[ERROR] Off BED could not be read : %s" % (args[1])
		exit(1)

	# read putative on-target region file (triplexator tts file)
	targets = {}
	try:
		for line in fileinput.input([args[2]]):
			if (len(line.strip())==0 or line[0]=="#"):
				continue
			cols = line.rstrip('\n').split("\t")
			if (not targets.has_key(cols[6])):
				targets[cols[6]] = []
			targets[cols[6]] += [cols]
		if (options.verbose):
			print >> sys.stderr, "TTS: read %d loci of interest"  % (len(targets))
	except:
		print >> sys.stderr, "[ERROR] primary target file could not be read : %s" % (args[2])
		exit(1)
		
	# read off-targets file (triplexator tpx file augmented with annotation data)
	tfotargets = {}
	signature = [] # the signature provides information about the annotations used
	annotations = []
	try:
		for line in fileinput.input([args[3]]):
			# keep signature
			if (line.startswith("#")):
				signature = line.strip("#").strip().split("\t")
				continue
	
			elif (len(line.strip())==0):
				continue
				
			cols = line.strip().split("\t")
			tfoId = cols[0]
	
			tfo_errors = ["-" for i in range(int(cols[1]),int(cols[2]))]
			offLength = int(cols[5])-int(cols[4])
			
			# requires the error to be with respect to the TFO
			errors = cols[8]
			if (errors!=""):
				i = 0
				while (i<len(errors)):
					if (errors[i] in "dobt"):
						# get number
						j=i+1
						while(j<len(errors) and errors[j] in "0123456789"):
							j += 1 
						tfo_errors[int(errors[i+1:j])] = errors[i]  
						i=j
			tfo_errors = ''.join(tfo_errors)
			
			# add new TFO entry if required
			if (not tfotargets.has_key(tfoId)):
				tfotargets[tfoId] = {}
				
			# add new tfo_region if required
			tfo_region = tuple([int(cols[1]),int(cols[2])])
			if (not tfotargets[tfoId].has_key(tfo_region)):
				tfotargets[tfoId][tfo_region] = {}
				
			# add new error category (string) if required
			if (not tfotargets[tfoId][tfo_region].has_key(tfo_errors)):
				tfotargets[tfoId][tfo_region][tfo_errors] = []
			# add relevant off-target information
			tfotargets[tfoId][tfo_region][tfo_errors] += [tuple([cols[3]]+[int(cols[4])]+[int(cols[5])]+cols[12:])]
	
		if (options.verbose):
			print >> sys.stderr, "TRIPLEX.file read (%d entries)"  % (len(tfotargets))
	except:
		print >> sys.stderr, "[ERROR] off-target file could not be read : %s" % (args[2])
		traceback.print_exc(file=sys.stderr)
		exit(1)

	#sort
	for tfoId in tfotargets.keys():
		for tfo_region in tfotargets[tfoId].keys():
			for tfo_errors in tfotargets[tfoId][tfo_region].keys():
				tfotargets[tfoId][tfo_region][tfo_errors] = sorted(tfotargets[tfoId][tfo_region][tfo_errors], key=itemgetter(0,1))

	# annotations are located after position 12
	annotations = signature[12:]
	
	# check if chromatin file has been supplied and test file
	if (options.chromatin != "" and options.chromatinFormat == "bam"):
		import pysam # lazy import of required module
		# try open the bam file 
		try:
			samfile = pysam.Samfile(options.chromatin, "rb" )
			samfile.close()
		except:
			print >> sys.stderr, "[ERROR] chromatin file (bam) could show_t be read : %s" % (options.chromatin)
			traceback.print_exc()
			exit(1)
	elif (options.chromatin != "" and options.chromatinFormat == "bigwig"):
		from bx.intervals.io import GenomicIntervalReader # lazy import of required module
		from bx.bbi.bigwig_file import BigWigFile  # lazy import of required module
		# try open the bigwig file
		try:
			bigwig = BigWigFile(open(options.chromatin))
		except:
			print >> sys.stderr, "[ERROR] chromatin file (bigwig) could not be read : %s" % (options.chromatin)	
    			traceback.print_exc()
			exit(1)
		
#	return (records, loiBed, offBed, targets, tfotargets, annotations, samfile, bigwig);
	
	
def processOfftarget(off_targets, targetid, tfo_region, tfo_error, lflank, rflank, tsize, tstrand, samfile, bigwig):
	global records, loiBed, offBed, targets, tfotargets, annotations

	counter = 1
	jsonErrors = len(tfo_error)-tfo_error.count('-')
	locFullStr = [" " for i in range(lflank+tsize+rflank)]
	
	if (tstrand=='+'):
		jsonOffset = lflank+tfo_region[0]
		for i in range(len(tfo_error)):
			locFullStr[lflank+tfo_region[0]+i] = tfo_error[i];
	else:
		jsonOffset = lflank+tsize-tfo_region[0]-len(tfo_error)
		for i in range(len(tfo_error)):
			locFullStr[lflank+tsize-tfo_region[0]-i-1] = tfo_error[i];
	
	for offtarget_entry in tfotargets[targetid][tfo_region][tfo_error]:
		chrString = offtarget_entry[0] 
		locChr = offBed[chrString][0]
		locStart = int(offBed[chrString][1]) + offtarget_entry[1]
		locStop = int(offBed[chrString][1]) + offtarget_entry[2]
			
		annostring = ""
		annoarray = [ 0 for i in range(len(annotations)) ]
		for i in range(len(annotations)):
			value = 0
			if (offtarget_entry[3+i] != "-"):
				value = len(offtarget_entry[3+i].split(","))
			annostring += annotations[i]+':{0:3d} '.format(value)
			annoarray[i] = value
	
		encoding = "".join(locFullStr).strip();
			
		off_chromatin_features = 0.
		if (options.chromatin != "" and options.chromatinFormat=="bam"):
			try:
				for pileupcolumn in samfile.pileup(locChr,locStart,locStop):
					off_chromatin_features += pileupcolumn.n
				off_chromatin_features /= (locStop-locStart)
			except:
				print >> sys.stderr, "[WARN] bam file exception : %s:%d-%d" % (locChr,locStart,locStop)
		elif (options.chromatin != "" and options.chromatinFormat=="bigwig"):
			try:
				bwsummary = bigwig.query(locChr, locStart, locStop, 1 ) 
				off_chromatin_features = bwsummary[0]["mean"]
			except:
				print >> sys.stderr, "[WARN] bigwig file exception : %s:%d-%d" % (locChr,locStart,locStop)
		else:
			off_chromatin_features = '-'

		if (not off_targets.has_key(jsonOffset)):
			off_targets[jsonOffset] = {}
			
		if (not off_targets[jsonOffset].has_key(encoding)):
			off_targets[jsonOffset][encoding] = {"errors":jsonErrors, "length":(locStop-locStart),"data":[],"copies":0, "chromatin_max":off_chromatin_features,"annotation":[0 for i in range(len(annotations))]}
			
		for anno in range(len(annotations)):
			off_targets[jsonOffset][encoding]["annotation"][anno] += annoarray[anno]

		if (off_chromatin_features > off_targets[jsonOffset][encoding]["chromatin_max"]):
			off_targets[jsonOffset][encoding]["chromatin_max"] = off_chromatin_features

		off_targets[jsonOffset][encoding]["copies"] += 1
		# skip display of off-targets > specified threshold
		if (len(off_targets[jsonOffset][encoding]["data"]) <= options.threshold):
#			if (off_chromatin_features == '-'):
#				off_targets[jsonOffset][encoding]["data"] += [(locChr, locStart, locStop, off_chromatin_features,'"'+'" "'.join(offtarget_entry[3:]).replace(',',' ').replace('" "','","')+'"')] 
#			else: # convert chromatin entry to string
#				off_targets[jsonOffset][encoding]["data"] += [(locChr, locStart, locStop, "%.3f" % (off_chromatin_features),'"'+'" "'.join(offtarget_entry[3:]).replace(',',' ').replace('" "','","')+'"')]
			anno_overlap = []
			for ai in offtarget_entry[3:]:
				astr = '-';
				if (ai != '-'):
					astr = str(len(ai.split(',')))
				anno_overlap += [astr]

			if (off_chromatin_features == '-'):
				off_targets[jsonOffset][encoding]["data"] += [(locChr, locStart, locStop, off_chromatin_features, '"'+'","'.join(anno_overlap)+'"' )]
			else: # convert chromatin entry to string
				off_targets[jsonOffset][encoding]["data"] += [(locChr, locStart, locStop, "%.3f" % (off_chromatin_features),'"'+'","'.join(anno_overlap)+'"' )]

		counter += 1 


def process():
	# get the data
	readdata()
	
	# set chromatin flag
	withChromatin = 'false'
	if (options.chromatin!=""):
		withChromatin = 'true'
		
	samfile = ""
	bigwig = ""
	# check if chromatin file has been supplied
	if (options.chromatin != "" and options.chromatinFormat == "bam"):
		import pysam # lazy import of required module
		# try open the bam file 
		try:
			samfile = pysam.Samfile(options.chromatin, "rb" )
		except:
			exit(1)
	elif (options.chromatin != "" and options.chromatinFormat == "bigwig"):
		from bx.intervals.io import GenomicIntervalReader # lazy import of required module
		from bx.bbi.bigwig_file import BigWigFile  # lazy import of required module
		# try open the bigwig file
		try:
			bigwig = BigWigFile(open(options.chromatin))
		except:
			exit(1)

	processors = 1
	if (options.processors <= 0): # figure out number of processors
		try:
			processors = multiprocessing.cpu_count()
		except:
			processors = 1
	else:
		processors = options.processors

	output_ontargets = "" # gets filled with a json object for on-targets
	
	# process one on-target region/locus of interest (LOI) at a time
	for loi in targets.keys():
		# info for whole region
		bed = loiBed[loi]
		chrom = bed[0]
		cstart = int(bed[1])
		cstop = int(bed[2])
		clabel = bed[3]
		cstrand = bed[5]
		tts = records[loi].seq	

		# process individual eligible target sites in this region/locus
		for target in targets[loi]:
			tstart = int(target[1])
			tstop = int(target[2])
			targetid = target[3]
			tstrand = target[5]
			tsize = tstop-tstart
			lflank = min(tstart-cstart, options.flanks)			
			rflank = min(cstop-tstop, options.flanks)
			
			targetSeq = ""
			if (cstrand == "+"):
				targetSeq = tts[tstart-cstart-lflank: tstop-cstart+rflank]
			else:
				targetSeq = tts.reverse_complement()[tstart-cstart-lflank: tstop-cstart+rflank]
			targetStr = [" " for i in range(lflank+tsize+rflank)]
			
			for i in range(lflank, lflank+tsize):
				targetStr[i]='+'
				
			if (target[8]!=""):
				for s in target[8].strip('d').split('d'):
					if (cstrand == "+"):
						targetStr[lflank+int(s)]='_'
					else:
						targetStr[lflank+tsize-int(s)-1]='_'
			
			on_chromatin_features = 0.
			if (options.chromatin != "" and options.chromatinFormat=="bam"):
				try:
					for pileupcolumn in samfile.pileup(chrom,tstart-lflank,tstop+rflank):
						on_chromatin_features += pileupcolumn.n
					on_chromatin_features = "%.3f" % (on_chromatin_features/(lflank+tsize+rflank))
				except:
					print >> sys.stderr, "[WARN] bam file exception : %s:%d-%d" % (chrom,tstart-lflank,tstop+rflank)
			elif (options.chromatin != "" and options.chromatinFormat=="bigwig"):
				try:
					bwsummary = bigwig.query(chrom, tstart-lflank, tstop+rflank, 1 ) 
					on_chromatin_features = "%.3f" % (bwsummary[0]["mean"])		
				except:
					print >> sys.stderr, "[WARN] bigwig file exception : %s:%d-%d" % (chrom,tstart-lflank,tstop+rflank)
			else:
				on_chromatin_features = "-"
			
			off_targets = {} # hashing start positions (w/r/t on target) and encoding (grouping off-targets) 				
			
			# branch point
			if (processors == 1): # no need to multiprocess

				for tfo_region in tfotargets[targetid].keys():
					# process each error category
					for tfo_error in tfotargets[targetid][tfo_region].keys():
						processOfftarget(off_targets, targetid, tfo_region, tfo_error, lflank, rflank, tsize, tstrand, samfile, bigwig);
						
			else: # contribute workload
				# load up work queue
				work_queue = multiprocessing.JoinableQueue(processors*3)
				
				# create a queue to pass to workers to store the results
				result_queue = multiprocessing.Queue()
				
				num_jobs = 0
				for tfo_region in tfotargets[targetid].keys():
					num_jobs += len(tfotargets[targetid][tfo_region])
					
				# spawn workers
				for i in range(min(num_jobs,processors)):
					worker = Worker(work_queue, result_queue, options.chromatin, options.chromatinFormat)
					worker.daemon=True
					worker.start()
					
				# process off-targets w/r/t to on-target
				# process each intersecting interval
				for tfo_region in tfotargets[targetid].keys():
					# process each error category
					for tfo_error in tfotargets[targetid][tfo_region].keys():
						work_queue.put(Job(targetid, tfo_region, tfo_error, lflank, rflank, tsize, tstrand));
						
				work_queue.close()		
				work_queue.join()
	
				for job in range(num_jobs):
					result = result_queue.get()
					for offset in result.keys():
						if (not off_targets.has_key(offset)):
							off_targets[offset] = {}
						for encoding in result[offset].keys():
							if (not off_targets[offset].has_key(encoding)):
								off_targets[offset][encoding] = result[offset][encoding]
			
			# merge point			
			output_sequence_w = ""
			output_sequence_c = ""
			for i in range(len(targetStr)):
				if (targetStr[i] ==" "):
					output_sequence_w += targetSeq.lower()[i]
					output_sequence_c += targetSeq.complement().lower()[i]
				else:
					output_sequence_w += targetSeq.upper()[i]
					output_sequence_c += targetSeq.complement().upper()[i]				

			# generate oaColumns (data table header)
			output_offtargets = ""
			for i in range(len(annotations)):
				output_offtargets += (',\n		{ "sTitle":"%s", "sToolTip":"number of off-targets intersecting with annotation data: %s" }' % (annotations[i],annotations[i]))
			output_offtargets += (',\n		{ "sTitle":"","bSortable":false,"bSearchable":false,"sClass":"btn_details","sToolTip":"Toggle detailed off-target inspector" }\n	],\n	"aaData":[')
			

			# generate aaData (data table content)
			for offset in off_targets.keys():
				for encoding in off_targets[offset]:
					ots = off_targets[offset][encoding]
					# format chromatin data if available
					if (options.chromatin != ""):
						ots["chromatin"] = "%.3f" % (ots["chromatin_max"])
					else:
						ots["chromatin"] = "-"

					output_offtargets += '\n		["%s",0,%d,%d,"%d-%d",%d,"%s",%d,[' % (encoding, ots["errors"], offset, offset-lflank, offset+len(encoding)-lflank,ots["copies"], ots["chromatin"], ots["length"],)
					for ot in ots["data"]:
						output_offtargets += '\n			["%s",%d,%d,"%s",%s],' % (ot)
					
					output_offtargets = output_offtargets[:-1]+"\n		]"
					if (len(annotations)>0):
						output_offtargets += ",%s" % (",".join(["%s" % el for el in ots["annotation"]]))
					output_offtargets += ',""],' 
			
			# finalise off-target json object and write output
			output_offtargets = ('''{
	"bPaginate": true,
	"bProcessing": true,
	"bAutoWidth": false,
	"bInfo": true,
	"bLengthChange": true,
    "bFilter": true,
	"iDisplayLength": 10,
	"aLengthMenu": [[10, 50, 100, -1], [10, 50, 100, "All"]],
	"bJQueryUI": true,
	"bDeferRender": true,
	"aaSorting": [[2, "asc"]],
	"aoColumns": [
		{ "sTitle":"%s","sClass":"monospace left", "sToolTip":"shows the part of the primary target that is responsible for off-targets" },
		{ "sTitle":"overlap", "sType": "numeric","sToolTip":"number of nucleotide positions off-targets overlap the primary target" },
		{ "sTitle":"errors", "sType": "numeric","sToolTip":"expected number of mismatches between the triplex-forming molecule designed against the primary target and off-targets" },
		{ "sTitle":"offset","sType": "numeric","bSearchable":false,"bVisible":false, "sToolTip":"the offset is used internally for computing the lefthand-side alignment"  },
		{ "sTitle":"sub region", "sToolTip":"region of the primary target this off-target spans" },
		{ "sTitle":"copies", "sType": "numeric","sToolTip":"number of copies of this off-target category" },
		{ "sTitle":"max. chromatin", "bVisible": %s, "sToolTip":"maximal chromatin score observed in any of the off-targets" },
		{ "sTitle":"length", "sType": "numeric","sToolTip":"length of the off-target", "bVisible":false },
		{ "sTitle":"offtargets","bVisible": false, "sToolTip":"detailed off-target list with location and annotation"  } '''+output_offtargets[:-1]+"\n	]\n}") % (output_sequence_w+"<br/>"+output_sequence_c, withChromatin)

		
			output = open(options.output_dir+targetid+"_off_targets.json","w")
			output.write(output_offtargets)
			output.close()	
			
			# pretty sequence encoding for on-target region
			output_sequence_w = ""
			output_sequence_c = ""
			lastpos = " "
			interruptions = 0
			for i in range(len(targetStr)):
				if (targetStr[i] != lastpos):
					if (lastpos !=" "):
						output_sequence_w+="</span>"
						output_sequence_c+="</span>"
					if (targetStr[i] =="+"):
						output_sequence_w += "<span class='valid_triad'>"
						output_sequence_c += "<span class='valid_triad'>"
					elif (targetStr[i] =="_"):
						output_sequence_w += "<span class='invalid_triad'>"
						output_sequence_c += "<span class='invalid_triad'>"
					lastpos = targetStr[i]
					if (targetStr[i] =="_"):
						interruptions += 1
				if (targetStr[i] ==" "):
					output_sequence_w += targetSeq.lower()[i]
					output_sequence_c += targetSeq.complement().lower()[i]
				else:
					output_sequence_w += targetSeq.upper()[i]
					output_sequence_c += targetSeq.complement().upper()[i]				
			if (lastpos == "+" or lastpos == "_"):
				output_sequence_w+="</span>"
				output_sequence_c+="</span>"
			output_sequence = "<span class='monospace flank_triad'>5'-"+output_sequence_w+"-3'<br/>3'-"+output_sequence_c+"-5'</span>"

			# count the number of off-targets
			offtarget_counter = 0
			for tfo_region in tfotargets[targetid].keys():
				for tfo_errors in tfotargets[targetid][tfo_region].keys():
					offtarget_counter+=len(tfotargets[targetid][tfo_region][tfo_errors])

			output_ontargets += '\n\t\t["%s","%s",%d,%d,%d,"%s",%d,"%s",%d,"%s",%d,%d],' % (targetid, chrom, tstart, tstop, (tstop-tstart), on_chromatin_features, interruptions, output_sequence, offtarget_counter, tstrand, lflank, rflank)
	
	# finalizes on-target json object and write output
	output_ontargets = '''{
	"bPaginate": true,
	"bProcessing": true,
	"bAutoWidth": false,
	"bInfo": true,
	"bLengthChange": true,
    "bFilter": false,
	"iDisplayLength": 10,
	"bJQueryUI": true,
	"aLengthMenu": [[5, 10, 25, -1], [5, 10, 25, "All"]],
	"aoColumns": [
		{ "sTitle": "region Id", "sToolTip":"internal id used to identify regions of overlapping putative primary targets" },
		{ "sTitle": "chr", "sToolTip":"chromosome this region is located in" },
		{ "sTitle": "start", "sType": "numeric", "sToolTip":"chromosomal start position of this region" },
		{ "sTitle": "end", "sType": "numeric", "sToolTip":"chromosomal end position of this region" },
		{ "sTitle": "length", "sType": "numeric", "sToolTip":"number of nucleotides spanned by this target region" },
		{ "sTitle": "chromatin", "bVisible": '''+withChromatin+''', "sToolTip":"chromatin score averaged over the shown region" },
		{ "sTitle": "Y-interruptions", "sType": "numeric", "sToolTip":"number of pyrimdine interruptions in the polypurine/polypyrimidine tract of the region" },
		{ "sTitle": "on-target region", "bSortable": false, "sClass": "left", "sToolTip":"sequence of the region (plus additional flanking positions if specified)"  },
		{ "sTitle": "off-targets", "bSearchable": false, "bVisible": false, "sToolTip":"total number of off-targets accumulated over all putative primary targets in this region" },
		{ "sTitle": "strand", "bSearchable": false, "bVisible": false, "sToolTip":"strand on which the purine tract is located" },
		{ "sTitle": "offset left", "sType": "numeric", "bSearchable": false, "bVisible": false, "sToolTip":"number of upstream flanking positions shown in the sequence" },
		{ "sTitle": "offset right", "sType": "numeric", "bSearchable": false, "bVisible": false, "sToolTip":"number of downstream flanking positions shown in the sequence" }
	],	
	"aaData": [''' + output_ontargets[:-1] + '\n	]\n}' 
	
	output = open(options.output_dir+"primary_target_regions.json","w")
	output.write(output_ontargets)
	output.close()	
	

	# close file handle
	if (options.chromatin != "" and options.chromatinFormat == "bam"):
		try:
			samfile.close()
		except:
			exit(1)
			
# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] LOI.bed LOC.bed TTSpool.file TRIPLEX.file 

Generates the json files for primary target regions, 
incorporating chromatin data on demand.
	'''
	
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-s", "--fasta-file", type="string", dest="fastaFile", default="", 
					help="location with fasta sequences (overwrites any specified column)")
	parser.add_option("-c", "--chromatin", type="string", dest="chromatin", default="", 
					help="location of a chromatin data")
	parser.add_option("-F", "--chromatin-format", type="string", dest="chromatinFormat", default="", 
					help="either bam (requires pysam) or bigwig (requires bx-python)")
	parser.add_option("-f", "--flanks", type="int", dest="flanks", default=5, 
					help="flanking sequence shown around primary target sites")
	parser.add_option("-o", "--output-dir", type="string", dest="output_dir", default="", 
					help="directory where the output files will be saved into")
	parser.add_option("-t", "--threshold", type="int", dest="threshold", default=25,
					help="maximum number of off-targets to display for a target")
	parser.add_option("-p", "--processors", type="int", dest="processors", default=1,
					help="number of processors to use; 0 = determine automatically; default 1")

	(options, args) = parser.parse_args()
	if (len(args) != 4):
		parser.print_help()
		parser.error("incorrect number of arguments")
	
	if (options.verbose):
		print >> sys.stderr, "collect Output from \n(LOI.bed) : %s\n(LOC.bed) : %s\n(TTSpool.file) : %s\n(TRIPLEX.file) : %s" % (args[0],args[1],args[2],args[3])
		if (options.chromatin != ""):
			print >> sys.stderr, "chromatin data : %s \nformat : %s\n" % (options.chromatin, options.chromatinFormat)			
	process()

	
if __name__ == "__main__":
	main()
	
