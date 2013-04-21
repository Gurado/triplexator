#!/usr/bin/python

import sys, os, glob, traceback
from optparse import OptionParser
import fileinput
from operator import itemgetter

global options
global args

# read all the required data
def readdata():

	# read LOI BED file
	lois = {}
	for line in fileinput.input([args[0]]):
		if (len(line.strip())==0 or line[0]=="#"):
			continue
		cols = line.strip().split("\t")
		lois[cols[3]] = {}
		# get strand infoq
		strand =""
		if (len(cols)>=6):
			strand = cols[5]

		# get exons
		exons = []
		if (len(cols)>=11):
			sizes = cols[10].split(",")
			starts = cols[11].split(",")
			for exon in xrange(int(cols[9])):
				exons += [tuple([int(cols[1])+int(starts[exon]), int(cols[1])+int(starts[exon])+int(sizes[exon])])]

		# save data
		lois[cols[3]] = {'chrom':cols[0], 'start':int(cols[1]), 'end':int(cols[2]), 'targets':[], 'strand':strand, 'exons':exons, 'gid': cols[4]}
		
	if (options.verbose):
		print >> sys.stderr, "LOI BED: read %d entries"  % (len(lois))
	
	# the primary target site bed file provides length information about the primary sites
	tts = {}
	for line in fileinput.input([args[1]]):
		cols = line.strip().split("\t")
		if (len(cols)<4):
			continue
		try:
			tts[cols[3]] = {'chrom':cols[0], 'start':int(cols[1]), 'end':int(cols[2])}
			source = cols[6]
			lois[source]['targets'] += [(cols[3], int(cols[1]))]
		except:
			if (options.verbose):
				print >> sys.stderr, 'skipping line TTS.file %s' % (line)

	# done
	if (options.verbose):
		print >> sys.stderr, "TTSs: read %d entries"  % (len(tts))

	
	# sort primary targets subject to start coordinate
	for loi in lois.keys():
		lois[loi]['targets'] = sorted(lois[loi]['targets'], key=itemgetter(0))

	# the tpx file provides all information about the off targets
	tpx = {}
	signature = [] # the signature provides information about the annotation used
	tpx_counter = 0
	for line in fileinput.input([args[2]]):
		cols = []
		# keep signature
		if (line.startswith("#")):
			signature = line.strip("#").strip().split("\t")
			continue
		
		cols = line.strip().split("\t")
		if (len(cols)<6):
			if (options.verbose):
				print >> sys.stderr, 'skipping line %s' % (line)
			continue
		
		ht = {}	
		for i in range(len(cols)):
			key = signature[i]
			value = ""
			try:
				value = int(cols[i])
			except:
				try:
					value = float(cols[i])
				except:
					value = cols[i]
			ht[key]=value
			
		oId = ht["Sequence-ID"]
		score = ht["Score"]

		# check if Sequence-ID is known
		if (not tts.has_key(oId)):
			print >> sys.stderr, 'Sequence-ID not found in TTS-file %s' % (line)
			continue
			
		# register entry
		if (not tpx.has_key(oId)):
			tpx[oId] = []

		tpx[oId] += [ht]
		tpx_counter += 1

	# done
	if (options.verbose):
		print >> sys.stderr, "TPXs: read %d entries"  % (tpx_counter)


	# process and save log files	
	parameters = {}
	tpx_logfile = args[2]+".log"
	try:
		if (options.tpx_logfile != ""):
			tpx_logfile = options.tpx_logfile
		log = open(tpx_logfile).read()
		parameters["errorRate"] = int(log.split("- maximum error-rate :")[1].split("%")[0])
		parameters["maxError"] = int(log.split("- maximum total error :")[1].split("\n")[0])
		parameters["minLength"] = int(log.split("- minimum length :")[1].split("nucleotides")[0])
		parameters["minGuanineRate"] = int(log.split("- minimum guanine content with respect to the target : ")[1].split("%")[0])
		parameters["maxGuanineRate"] = int(log.split("- maximum guanine content with respect to the target : ")[1].split("%")[0])
		parameters["matchBlock"] = int(log.split("- number of consecutive matches required in a feature :")[1].split("\n")[0])
	except:
		print >> sys.stderr, 'log file not found or in unexpected format: %s' % (tpx_logfile)
		traceback.print_exc()
		sys.exit(1)
	options.tpx_logfile = os.path.abspath(tpx_logfile)
	# make relative to output dir
	if (options.output != ""):
		options.tpx_logfile = options.tpx_logfile.split(os.path.abspath(options.output))[-1].lstrip('/')

	try:
		tts_logfile = args[1].rstrip(".bed").rstrip("pool")+'.log'
		if (options.tts_logfile == ""):
			options.tts_logfile = os.path.abspath(tts_logfile)
		else:
			options.tts_logfile = os.path.abspath(options.tts_logfile)
		# make relative to output dir
		if (options.output != ""):
			options.tts_logfile = options.tts_logfile.split(os.path.abspath(options.output))[-1].lstrip('/')
	except:
		print >> sys.stderr, 'primary target log file not found: %s' % (options.tts_logfile)
		sys.exit(1)
		
	try:
		if (options.workflow_logfile != "log.txt"):
			options.workflow_logfile = os.path.abspath(options.workflow_logfile)
			# make relative to output dir
			if (options.output != ""):
				options.workflow_logfile = options.workflow_logfile.split(os.path.abspath(options.output))[-1].lstrip('/')
	except:
		print >> sys.stderr, 'workflow log file not found: %s' % (options.workflow_logfile)
		sys.exit(1)
	
	# done
	if (options.verbose):
		print >> sys.stderr, "...finished"

	return (signature, tpx, tts, lois, parameters)

def write_html(signature, tpx, tts, lois, parameters):

	html_header = '''<!DOCTYPE html>
	<html>
	<head>
	<title>Triplex-Inspector</title>
	<meta name="author" content="Fabian Buske, The University of Queensland, Australia" />
	<link rel="shortcut icon" href="includes/image/favicon.ico" />
	<link rel="icon" href="includes/image/favicon.ico" type="image/x-icon" />
	
	<!-- include local version of jquery to support offline functionality -->
	<script type="text/javascript" src="includes/js/jquery.min.js"></script>
	<script type="text/javascript" src="includes/js/jquery-ui.min.js"></script>

	<!-- data tables -->
	<script type="text/javascript" src="includes/js/jquery.dataTables.min.js"></script>
	<script type="text/javascript" src="includes/js/jquery.dataTables.selectFilter.js"></script>
	<script type="text/javascript" src="includes/js/ColVis.min.js"></script>
	<script type="text/javascript" src="includes/js/FixedHeader.min.js"></script>
	<script type="text/javascript" src="includes/js/CustomHashTable.js"></script>

	<!-- tooltip -->
	<script type="text/javascript" src="includes/js/jquery.tools.min.js"></script>

	<!-- flotr2 -->
	<script type="text/javascript" src="includes/js/flotr2.min.js"></script>	

	<!-- style sheets -->
	<style type="text/css" title="currentStyle">
		@import "includes/css/inspector.css";
		@import "includes/css/jquery-ui.css";
		@import "includes/css/jquery.dataTables.css";
		@import "includes/css/ColVis.css";
		@import "includes/css/flotr2.css";
	</style>
	
	<!-- triplex inspector specific scripts -->
	<script type="text/javascript">
		var hasChromatin = %d;
		var hasUCSCgenome = %d;
		var UCSCgenome = "%s";
		var submatches_file = "%s";
		var parameters = {"minLength": %d, "errorRate": %d, "maxError": %d, "minGuanineRate": %d, "maxGuanineRate": %d, "matchBlock": %d}
    </script>
   	<script type="text/javascript" src="includes/js/inspector.js"></script>''' % (options.chromatin!="NONE", options.ucsc!="NONE", options.ucsc, args[3], parameters["minLength"], parameters["errorRate"], parameters["maxError"], parameters["minGuanineRate"], parameters["maxGuanineRate"], parameters["matchBlock"])
	
	report_body = '''
	<script type="text/javascript" src="includes/js/inspector_report.js"></script>
	</head> 
	<body>
	<div class="wrapper">

	<div id="target_region_head">Targetable regions</div>
	<div id="target_region_overview">
		<table id="tbl_sequences">
		<colgroup>
			<col style="width:20em;" />
			<col />
		</colgroup>
		<thead>
			<tr>
				<th>Genomic locus</th>
				<th>Target Region Block Diagram</th>
			</tr>
		</thead>
		<tbody>
	'''
	
	locationArray = "var locations = new Array();\n";
	
	for loi in lois.keys():
		loi_entry = lois[loi]
		loi_chrom = loi_entry['chrom']
		loi_start = loi_entry['start']
		loi_end = loi_entry['end']
		loi_strand = loi_entry['strand']
		loi_exons = loi_entry['exons']
		loi_size = float(loi_end - loi_start)
		loi_gid = loi_entry['gid']
		
		# calculate maximum stacking
		step = 1
		maxstep = step
		pre_tts_chrom = ""
		pre_tts_start = -1
		pre_tts_end = -1
		for target in loi_entry['targets']:
			primary_target = target[0] 
			tts_chrom = tts[primary_target]['chrom']
			tts_start = tts[primary_target]['start']
			tts_end = tts[primary_target]['end']
			if (pre_tts_chrom == tts_chrom and (pre_tts_start <= tts_start and tts_start <= pre_tts_end) or (pre_tts_start <= tts_end and tts_end <= pre_tts_end)):
				step += 1 
			else:
				step = 1
			pre_tts_chrom = tts_chrom
			pre_tts_start = tts_start
			pre_tts_end = tts_end
			maxstep = max(maxstep, step)
		
		if (options.ucsc!="NONE"):
			# print sequence block container
			report_body += '''
			<tr id='%s_blocks'>
				<td><a href='http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s:%d-%d' target='ucsc'>%s:%d-%d (%s)</a><br/><span class="small">%s</span></td>
				<td>
					<div id='%s_block_container' class='block_container' style='height:%dpx'>
					<div class='block_rule' style='width:100%%'></div>
			''' % (loi, options.ucsc, loi_chrom, loi_start, loi_end, loi_chrom, loi_start, loi_end, loi_strand, loi, loi, maxstep*13+13)
		
		else:
			# print sequence block container
			report_body += '''
			<tr id='%s_blocks'>
				<td>%s:%d-%d (%s)<br/><span class="small">%s</span></td>
				<td>
					<div id='%s_block_container' class='block_container' style='height:%dpx'>
					<div class='block_rule' style='width:100%%'></div>
			''' % (loi, loi_chrom, loi_start, loi_end, loi_strand, loi, loi, maxstep*13+13)
	
		# print exon blocks
		for exon_i in xrange(len(loi_exons)):
			exon = loi_exons[exon_i]
			exon_nm = "e_"+str(exon_i)
			if (loi_strand == '-'):
				exon_nm = "e_"+str(len(loi_exons)-exon_i)
			report_body += '''		<div class='block_exon' id="%s_exon" title='%s %s:%d-%d' style='left:%.2f%%; top:%dpx; width:%.2f%%; border: 1px dotted black;'></div>''' % (exon_nm, exon_nm, loi_chrom, exon[0], exon[1], (exon[0]-loi_start)*100./loi_size, 8*step, (exon[1]-exon[0])*100./loi_size+0.1)  # block
		
		# print motif blocks
		step = 1
		pre_tts_chrom = ""
		pre_tts_start = -1
		pre_tts_end = -1
		for target in loi_entry['targets']:
			primary_target = target[0] 
			tts_chrom = tts[primary_target]['chrom']
			tts_start = tts[primary_target]['start']
			tts_end = tts[primary_target]['end']
			
			if (pre_tts_chrom == tts_chrom and (pre_tts_start <= tts_start and tts_start <= pre_tts_end) or (pre_tts_start <= tts_end and tts_end <= pre_tts_end)):
				step += 1 
			else:
				step = 1

			pre_tts_chrom = tts_chrom
			pre_tts_start = tts_start
			pre_tts_end = tts_end
			
			report_body += '''		<div class='block_motif' id="%s_block" title='%s %s:%d-%d' style='left:%.2f%%; top:%dpx; width:%.2f%%; border: 1px solid black;' onClick='showDetail("%s");'></div>''' % (primary_target, primary_target, tts_chrom, tts_start, tts_end, (tts_start-loi_start)*100./loi_size, 7*step, (tts_end-tts_start)*100./loi_size+0.1, primary_target)  # block
			report_body += '''		<div class='block_id' style='left:%.2f%%; top:%dpx; width:%.2f%%;'>%s</div>''' % ((tts_start-loi_start)*100./loi_size, 17*step, (tts_end-tts_start)*100./loi_size+0.1, primary_target)  # block

			
			locationArray += "	locations['%s'] = '%s:%d-%d';\n" % (primary_target, tts_chrom, tts_start, tts_end)
			
		# close table
		report_body += '''
			</div></td>
			</tr>'''
		
	report_body += '''
		</tbody>
	</table>
	
	<table cellpadding="0" cellspacing="0" border="0" class="display" id="onregion_table">
	<thead>
	<tr>
	</tr>
	</thead>
	<tbody>
	<tr>
	</tr>
	</tbody>
	</table>
	</div>

	<div id="on_target_head">Pool of putative targets</div>
	<div id="on_target_bowser">
	<table class="display" id="filter_table">
	    <thead>
	    	<tr>
	    		<th colspan="8" style="background-color:#DDFFDD">primary targets &amp; properties</th>
	    		<th colspan="5" style="background-color:#FFDDDD">off-target abundance</th>
	    	</tr>
	        <tr>
	            <th></th>
	            <th></th>
	            <th></th>
	            <th></th>
	            <th></th>
	            <th></th>
	            <th></th>
	            <th></th>
	            <th></th>
	            <th></th>
	            <th></th>
	            <th></th>
	            <th></th>
	        </tr>
	    </thead>
 		<tfoot>
	        <tr>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	            <th style="vertical-align: top; padding-top: 10px"></th>
	        </tr>
	    </tfoot>
	    <tbody>
	        <tr>
	            <td></td>
	            <td></td>
	            <td></td>
	            <td></td>
	            <td></td>
	            <td></td>
	            <td></td>
	            <td></td>
	            <td></td>
	            <td></td>
	            <td></td>
	            <td></td>
	            <td></td>
	        </tr>
	    </tbody>
	</table>
	<table border="0" class="small"><tr><td>Color-coding indicates the number of &#39;as good&#39; off-targets: </td>
	<td style="background-color:#ddffdd;width:15px;">&nbsp;</td><td>0-10</td>
	<td style="background-color:#ddddff;width:15px;">&nbsp;</td><td>11-100</td>
	<td style="background-color:#ffdddd;width:15px;">&nbsp;</td><td>101-1000</td>
	<td style="background-color:#dddddd;width:15px;">&nbsp;</td><td>&gt;1000</td></tr></table>
	</div>'''
	
	report_body += '''
	<div id="log_head">Log files</div>
	<div id="log_browser">
	<table>
	<tr><td>workflow log:</td><td><a href="%s">%s</a></td></tr>
	<tr><td>triplexator (primary target)</td><td><a href="%s">%s</a></td></tr>
	<tr><td>triplexator (off-targets)</td><td><a href="%s">%s</a></td></tr>
	</table>
	</div>
	''' % (options.workflow_logfile, options.workflow_logfile, options.tts_logfile, options.tts_logfile, options.tpx_logfile, options.tpx_logfile)
	
	detail_body = '''
	<script type="text/javascript" src="includes/js/inspector_detail.js"></script>
	</head> 
	<body>
	<div class="wrapper">

	<div id="off_target_head">Off-targets for primary target at <span id="target_details">-</span> <div style="float:right"><a href="inspector_report.html" style="color:#fff">Back to report</a></div></div>
	<div id="off_target_inspector">
	<table class="display" id="offtarget_table">
	    <thead>
	    	<tr>
	        </tr>
	    </thead>
	    <tfoot>
	    	<tr>
	        </tr>
	    </tfoot>
	    <tbody>
	        <tr>
	        </tr>
	    </tbody>
	</table>

	<!-- Mismatch/Risk profile + TFO sequences -->
	<table>
	<tr><td width="450">
		<div id="mismatches" style="width: 450px; height: 270px; margin: 8px;-moz-user-select: none;"></div>
	</td><td width="10">&nbsp;</td><td width="550">
		<div id="tfo" style="margin: 8px; visibility:hidden">
		<table border="0" id="tfo-b" width="100%">
		<thead>
			<tr>
				<th>Motif<sup>1</sup></th>
				<th>Proposed TFO sequence<sup>2</sup> for this target site</th>
				<th>Preferred<sup>3</sup></th>
			</tr>
		</thead>
		<tbody>
			<tr> 
				<td>TM</td>
				<td><span id="TMmotif" class="monospace"></span></td>
				<td align="center"><span id="TMmotif_preferred"></span></td>
			</tr>
			<tr>
				<td>GU</td>
				<td><span id="GUmotif" class="monospace"></span></td>
				<td align="center"><span id="GUmotif_preferred"></span></td>				
			</tr>
			<tr>
				<td>GA</td>
				<td><span id="GAmotif" class="monospace"></span></td>
				<td align="center"><span id="GAmotif_preferred"></span></td>
			</tr>
		</tbody>
		<tfoot>
			<tr>
				<td colspan="3" class="legend">
				<sup>1</sup> : TM-motif TFO containing thymidine (T) and 5-methyldeoxycytidine (M) - parallel binding<br/>	
				&nbsp;&nbsp;&nbsp;&nbsp; GU-motif TFO containing deoxyguanosine (G) and deoxyuridine (U) - anti-parallel binding<br/>
				&nbsp;&nbsp;&nbsp;&nbsp; GA-motif TFO containing deoxyguanosine (G) and deoxyadenosine (A) - anti-parallel binding<br/>
				<sup>2</sup> : TFO sequence templates according to the models by <a href="http://pubs.acs.org/doi/abs/10.1021/bi801087g" target="blank">Vekhoof <i>et al.</i><a/> <br/>
				<sup>3</sup> : Preference is calculated using formula (1) in <a href="http://pubs.acs.org/doi/abs/10.1021/bi801087g" target="blank">Vekhoof <i>et al.</i></a> and is based on assumptions given therein.<br/>
				<span class="rank0">x</span> : position requires nucleotide choice due to a pyrimidine interruption in the primary target.<br/>
				<span class="rank1">n</span> : position a strategically placed mismatch will affect the most off-target sites (absolute).<br/>

				</td>
			</tr>
		</tfoot>
		</table>
		</div>
	</td></tr></table>

	<u>Legend:</u> with respect to the <a href='http://acb.qfab.org/acb/triplexator/biology.html' target='_blank'>canonical triplex formation ruleset</a><br/>
	'-': valid triplex triad formed between the third strand and the target duplex at this position (primary and off-target respectively), <br/>
	'o': error due to the triplex-forming molecule (third strand) as consequence of an error in the primary target, <br/> 
	'd': error due to the off-target duplex (pyrimidine interruption), <br/>
	'b': error due to both the triplex-forming molecule and the duplex,<br/>
	't': error due to the triplex pairing between the third strand nucleotide and the duplex,<br/>
	'x': pyrimidine interruption in the primary target		
	</div>
	'''
	
	html_footer = '''<div id='footer'>Copyright &copy; 2012, <a href='mailto:fabian.buske@gmail.com'>Fabian Buske</a>. All rights reserved.</a> Check <a href='./FAQs.txt'>FAQs</a> when encountering issues.</div></div>
	<div id="browsererror" class="center" style="position:absolute;top:100px;"><div class="error">If you can read this, check the <a href="FAQs.txt">FAQs</a>.</div></div>
	<div id="datawarning" class="center" style="position:absolute;top:100px;visibility:hidden"><div class="warning">No region has been specified. Go to the report and select a target.</div></div>
	<div id="datanotice" class="center" style="position:absolute;top:100px;visibility:hidden"><div class="notice">Obtain data from server. Please wait</div></div>

</body>\n</html>'''

	report = open(options.output+"inspector_report.html","w");
	report.write(html_header)
	report.write(report_body)
	report.write(html_footer)
	report.close()
	detail = open(options.output+"inspector_detail.html","w");
	detail.write(html_header)
	detail.write(detail_body)
	detail.write(html_footer)
	detail.close()

def process():
	"""
	workflow
	"""
	
	(signature, tpx, tts, lois, parameters) = readdata()
	write_html(signature, tpx, tts, lois, parameters)
	
# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] loi.bed tts.bed tpx.file primaryTargets.json

takes two bed files containing 
(1) generell regions of interest, 
(2) location of the primary target region 
(3) triplexator output file containing information about the off-targets 
(4) location of the primary target json file
and generates a report in html format. 
	'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-o", "--output", type="string", dest="output", default="", 
					help="where to write the output files")
	parser.add_option("-x", "--off-target-log", type="string", dest="tpx_logfile", default="", 
					help="triplexator log file to retrieve additional parameters (default: looks for tpx.log file in the TFP.file location)")
	parser.add_option("-p", "--primary-target-log", type="string", dest="tts_logfile", default="", 
					help="triplexator log file for primary target detection (default: derived from tts.bed)")
	parser.add_option("-w", "--workflow-log", type="string", dest="workflow_logfile", default="log.txt", 
					help="location of the workflow log file w/r/t the output (default: log.txt)")
	parser.add_option("-c", "--with-chromatin", type="string", dest="chromatin", default="NONE", 
					help="whether to incorporate chromatin data")
	parser.add_option("-u", "--with-ucsc", type="string", dest="ucsc", default="NONE", 
					help="genome assembly to generate links to the UCSC genome browser")
	
	(options, args) = parser.parse_args()
	if (len(args) != 4):
		parser.print_help()
		parser.error("incorrect number of arguments")
		
	if (options.output != ""): 
		options.output += "/"
		

	if (options.verbose):
		print >> sys.stderr, "LOI file  : %s" % (args[0])
		print >> sys.stderr, "TTS file  : %s" % (args[1])
		print >> sys.stderr, "TPX file  : %s" % (args[2])
		print >> sys.stderr, "Json file  : %s" % (args[3])
			
	process()

	
if __name__ == "__main__":
	main()
