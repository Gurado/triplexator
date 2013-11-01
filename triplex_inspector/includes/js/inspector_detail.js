var dTable = null;
var onTargetRegion = null;
var offTargetRegionIds = [];
var annotations = [];
var chromatin = false;
var thisUrlVars;
var onOffsetLeft, onOffsetRight;
var onStart;
var onEnd;
var onSize;
var onStrand;
var mismatchdata = [];
var guaninePositions = []; // positions where guanines are in the primary target
var numberGuanines = 0; // number of guanines in the primary target
var purineTrack = ""; // target sequence containing the purine track
var TCmotif = "", GAmotif="", GTmotif="";
var maxDisruptingPos = new Array();


function bars_stacked() {
	Flotr.draw($('#mismatches')[0], mismatchdata, {
		legend: {
			position: "te"},
        bars: {
            show: true,
            stacked: true,
            horizontal: false,
            barWidth: 1,
            lineWidth: 1,
            shadowSize: 0
        },
		mouse: {
			track: true,
			relative: true
		},
		xaxis: {
			min: onOffsetLeft-0.5,
			max: onSize-onOffsetRight+0.2*onSize,
			title: "nucleotide position"
		},
		yaxis: {
			title: "risk (off-targets)",
		},
		title: "Positional off-target risk given parameter setting",
		spreadsheet: {show: true},
		HtmlText: false,
	});
}

/* calculate TFOs for the primary target using the canonical rules */
function tfo_motifs() {
	var tc = [],
		ga = [],
		gt = [];
	var leftp, rightp = 0;
	if (onStrand == '+'){
		leftp = onStart;
		rightp = onEnd;
	} else { // adjust pointers for crick strand targets
		leftp = onSize - onEnd;
		rightp = onSize - onStart;
	}
//	alert(leftp+' ' +rightp)
	for (var i=leftp; i<rightp; i++){
		if ($.inArray(i, maxDisruptingPos) != -1){
			leftMark = '<span class="rank1">';
			rightMark = '</span>';
		} else {
			leftMark = '';
			rightMark = '';		
		}
		switch(purineTrack[i]){
			case 'A':
				tc.push(leftMark+'T'+rightMark);
				ga.push(leftMark+'A'+rightMark);
				gt.push(leftMark+'U'+rightMark);
				break;
			case 'G':
				tc.push(leftMark+'M'+rightMark);
				ga.push(leftMark+'G'+rightMark);
				gt.push(leftMark+'G'+rightMark);
				break;
			default: // error in primary target
				tc.push('<span class="rank0">x</span>');
				ga.push('<span class="rank0">x</span>');
				gt.push('<span class="rank0">x</span>');							
		}
	}
	$('#TMmotif').html('5&prime;-'+tc.join('')+'-3&prime;');
	$('#GAmotif').html('5&prime;-'+ga.reverse().join('')+'-3&prime;');
	$('#GUmotif').html('5&prime;-'+gt.reverse().join('')+'-3&prime;');	
	$('#tfo').css('visibility','visible');
	
	var Gn = 73.0-1.1*(onEnd-onStart); //according to Vekhoff et al. 
	var Gr = 100*numberGuanines/(onEnd-onStart);
	if (Gr <= Gn){
		$('#TMmotif_preferred').html('&#10004;');
	} else {
		$('#GUmotif_preferred').html('&#10004;');
	}
}

/* Formating function for row details */
function fnFormatDetails ( dTable, nTr ){
    var aData = dTable.fnGetData( nTr );
	var sOut = '<table cellpadding="0" cellspacing="0" border="0" style="width:99%;">';
	sOut +=	'<thead><tr><th>off-target location</th><th>chromatin</th>'
	for (i=0;i<annotations.length;++i){
		sOut +=	'<th>'+annotations[i]+'</th>';
	}
	sOut +=	'</tr></thead>'

	for (i=0;i<aData[8].length;++i){
		if (hasUCSCgenome){
			sOut +=	'<tr><td align="left"><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+UCSCgenome+'&position='+aData[8][i][0]+':'+aData[8][i][1]+'-'+aData[8][i][2]+'" target="ucsc">'+aData[8][i][0]+':'+aData[8][i][1]+'-'+aData[8][i][2]+'</a></td>';
				for (j=3;j<aData[8][i].length;++j){
					sOut +=	'<td>'+aData[8][i][j]+'</td>';
				}
		} else {
			sOut +=	'<tr><td align="left">'+aData[8][i][0]+':'+aData[8][i][1]+'-'+aData[8][i][2]+'</td>';
				for (j=3;j<aData[8][i].length;++j){
					sOut +=	'<td>'+aData[8][i][j]+'</td>';
				}
		}
		sOut +=	'</tr>';
	}
	var copies = aData[5]
	if (copies > aData[8].length) {
		sOut +='<tr><td align="left" colspan="'+ (annotations.length+2)+'"><i>... +'+(copies-aData[8].length)+' more</i></td></tr>';	
	}
	sOut +='</table>';
    return sOut;
}

/* submatches with the same off-target signature are collapsed */
function fnCollapseOfftargetCategories(dTable){
	var rows = dTable.fnGetData();
	var leftString = new Array(onStart+1).join(' ');
	var signatures = new CustomHashTable(); // hold all different signatures
	var originalRowCount = rows.length;
	for(var i=0;i<originalRowCount;i++){
		signature = (leftString+$(rows[i])[0].substring(onStart,onEnd));
		if (! signatures.hasItem(signature)){
			signatures.setItem(signature, new Array());
		} 
		svalue = signatures.getItem(signature);
		svalue.push(i); // add row
		signatures.setItem(signature, svalue);
	}
	
	// add new elements
   	var entries = new Array();
	for (var k in signatures.items) {
	    if (signatures.hasItem(k)) {
	    	var entry = false;
//	    	alert('key is: ' + k + ', value is: ' + signatures.items[k]);
	    	var collapseRows = signatures.items[k];
	        for (var r =0; r<collapseRows.length;++r){
	        	var row = collapseRows[r];
	        	if (!entry){
	        		entry = rows[row] // take row as template
	        		entry[0] = k // modify signature
	        		entry[1] = k.trim().length; // modify overlap
	        		entry[2] = entry[1]-(k.trim().split(/-/g).length - 1); // modify errors
	        		// keep offset (entry[3]) 
					entry[4] = Math.max(parseInt(entry[4].split('-')[0]), onStart-onOffsetLeft)+'-'+Math.min(parseInt(entry[4].split('-')[1]), onEnd-onOffsetLeft);
					// keep copies (entry[5])
					// keep chromatin (entry[6])
					entry[7] = k.trim().length; // length is overlapping portion
					// keep specific sites (entry[8])
					// keep all annotations (entry[9+])
	        	} else {
	        		var collapse = rows[row];
	        		entry[5] += collapse[5]; // add copies
	        		// add all sites
	        		for (var j=0; j< collapse[8].length; ++j){
	        			entry[8].push(collapse[8][j]);
	        		}
	        		// add all annotations
	        		for (var i=9; i<collapse.length-1;++i){ // last element is button
	        			entry[i] += collapse[i];
	        		}
	        	}
	        }
	        if (entry){
        		entries.push(entry);
	        }
	    }
	}
	// clear table
	dTable.fnClearTable(false);
	
	// add new collapsed entries
	// check for suitable intersection
	var kept_entries = new Array();
	for (var entryindex=0; entryindex<entries.length;entryindex++){
		var aData=entries[entryindex];
	   	var overlap = aData[1];
	   	var qualifies = true;
	   	if (overlap < parameters["minLength"]){
	   		qualifies = false;
	    } else {
			var offStart = aData[3];
	    	var offEnd = offStart+aData[7];
	    	var from = onStart > offStart ? onStart : offStart;
	    	var to = onEnd < offEnd ? onEnd : offEnd;
	    	var errors = 0;
	    	var matchBlock = 0, tmp_mb = 0;
	    	var guanines = 0;
	    	for (i=from; i< to ;++i){
	    		if (aData[0][i] != '-'){
	    			++errors;
	    			if (tmp_mb>matchBlock){
	    				matchBlock = tmp_mb;
	    			}
	    			tmp_mb = 0;
	    		} else {
	    			++tmp_mb;	
	    		}
	    		guanines += guaninePositions[i];
	    	}
	    	if (tmp_mb>matchBlock){
	  				matchBlock = tmp_mb;
	  			}
	  			if (errors > parameters["maxError"] || errors*100./overlap > parameters["errorRate"] || matchBlock < parameters["matchBlock"] || guanines*100./overlap < parameters["minGuanineRate"] || guanines*100./overlap > parameters["maxGuanineRate"]){
	   			qualifies = false;
	    	}
	    }
	    
	    // add qualifying entries
	    if (qualifies){
	    	kept_entries.push(aData);
	    }
	}	
  	dTable.fnAddData(kept_entries, true);
}

/* submatches with the same off-target signature are collapsed */
$.fn.dataTableExt.oApi.fnCalculateMismatches = function (oSettings){
	mismatchdata = [],
	ismatch = [];
	var errorTFO = [],
		errorTTS = [],
		errorBoth = [],
		errorTriplex = [],
		disruptingPos = [];
	
	for(var i=0;i<onEnd;i++){
		ismatch.push([i,0]);
		errorTFO.push([i,0]);
		errorTTS.push([i,0]);
		errorBoth.push([i,0]);
		errorTriplex.push([i,0]);
		disruptingPos.push([i,0]);
	}
	for ( var row=0, iLen=oSettings.aiDisplay.length ; row<iLen ; row++ ){
  		var entry  = oSettings.aoData[ oSettings.aiDisplay[row] ]._aData
		var offStart = entry[3];
    	var offEnd = offStart+entry[7];
    	var from = onStart > offStart ? onStart : offStart;
    	var to = onEnd < offEnd ? onEnd : offEnd;
    	for (var i=from; i< to ;++i){
   			switch(entry[0][i]){
			case '-':
				ismatch[i][1] += entry[5]; // add copies
	   			disruptingPos[i][1] += entry[5]; // add copies
	   			break;
			case 'o':
	   			errorTFO[i][1] += entry[5]; // add copies
	  			break;
			case 'd':
				errorTTS[i][1] += entry[5]; // add copies
	   			break;
			case 'b':
				errorBoth[i][1] += entry[5]; // add copies
	   			break;
			case 't':
				errorTriplex[i][1] += entry[5]; // add copies
	   			break;
			}  //  default {}uncategorised
	   	}
	}
	
	mismatchdata = [{
		data:ismatch.slice(onStart,onEnd),
		label:"Match (-)"
	},{
		data:errorTFO.slice(onStart,onEnd),
		label:"Primary (o)"
	},{
		data:errorTTS.slice(onStart,onEnd),
		label:"Off-target (d)"
	},{
		data:errorBoth.slice(onStart,onEnd),
		label:"Both (b)"
	},{
		data:errorTriplex.slice(onStart,onEnd),
		label:"Triplex (t)"
	}];

	var disruptScore = 0;
	for (var i=onStart; i<onEnd; i++){
		if (disruptingPos[i][1] > disruptScore){
			disruptScore = disruptingPos[i][1];
			maxDisruptingPos = new Array();
			if (onStrand == '+'){
				maxDisruptingPos.push(i);
			} else {
				maxDisruptingPos.push(onSize-1-i);				
			}
		} else if (disruptingPos[i][1] == disruptScore){
			if (onStrand == '+'){
				maxDisruptingPos.push(i);
			} else {
				maxDisruptingPos.push(onSize-1-i);				
			}
		}
	}
	
	bars_stacked();
	tfo_motifs();
}

function handleAjaxError( xhr, textStatus, error ) {
    if ( textStatus === 'timeout' ) {
        alert( 'The server took too long to send the data.' );
    } else {
        alert( 'An error occurred on the server. Please try again in a minute.' );
    }
    dTable.fnProcessingIndicator( false );
}


// populate table when the document is ready
$(document).ready(function(){
	thisUrlVars = $.getUrlVars();
	
	// check if a region Id has been provided
	if ($.getUrlVar('rId') == undefined){
		$('#browsererror').css('visibility','hidden');
		$('#datawarning').css({opacity: 0.0, visibility: "visible"}).animate({opacity: 1.0});
	}
	else {
	
		// get ontarget entry
		$.ajax({
			url: 'json/primary_target_regions.json',
			async: false,
			dataType: 'json',
			success: function (json) {
				// json works, hide warning
				$('#browsererror').css('visibility','hidden');
				$('#datanotice').css({opacity: 0.0, visibility: "visible"}).animate({opacity: 1.0});
				for (i=0;i<json["aaData"].length;i++){
					if (json["aaData"][i][0] === $.getUrlVar('rId')){
						onTargetRegion = jQuery.extend({},json["aaData"][i]);
						break;
					}
				}
			}
	
		});
		setTimeout("loadOfftargets()",250);
	}
});

function loadOfftargets(){
		
	var tid_details = onTargetRegion[1] +':'+ onTargetRegion[2] +'-'+ onTargetRegion[3] + ' ('+ $.getUrlVar('rId')+')';
	
	if ($.getUrlVar('region') != ""){
		tid_details += "; subregion "+$.getUrlVar('region');
	}
	if ($.getUrlVar('annotation') != "-" && $.getUrlVar('annotation') != ""){
		tid_details += "; annotation "+$.getUrlVar('annotation');
	}
	$("#target_details").text(tid_details);
	if (hasUCSCgenome){
		$("#target_details").append(' (<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+UCSCgenome+'&position='+onTargetRegion[1]+':'+onTargetRegion[2]+'-'+onTargetRegion[3]+'" target="ucsc" style="color:#fff">UCSC</a>)');
	}
			
	// submatches_file needs to be set prior to execution
	$.ajax({
		url: 'json/'+$.getUrlVar('rId')+'_off_targets.json',
		async: true,
		dataType: 'json',
		success: function (json) {
			// get primary target specifics from table header
			onOffsetLeft = onTargetRegion[10];
			onOffsetRight = onTargetRegion[11];
			onStart = onOffsetLeft;
			onEnd = onOffsetLeft+onTargetRegion[4];
			onSize = onTargetRegion[4]+onOffsetLeft+onOffsetRight;// size of the ontarget region (+ flanks)
			onStrand = onTargetRegion[9];
			if ($.getUrlVar('region')!=""){
				onStart = parseInt($.getUrlVar('region').split("-")[0])+onOffsetLeft; // start + flank offset
   				onEnd = parseInt($.getUrlVar('region').split("-")[1])+onOffsetLeft;
			}

			// get hash map for column titles // Attention, problems with reserved words such as length are a given
			for (i=0; i < json["aoColumns"].length; ++i){
				hash = json["aoColumns"][i]["sTitle"];
				offTargetRegionIds[hash] = i;
				if (i>=9 && i < json["aoColumns"].length-1){
					annotations.push(hash);
				}
			}
			
			// array with guanine positions (G) indicated by 1 (0 otherwise)
			var duplex = json["aoColumns"][0]["sTitle"].split('<br/>')
			var strand = (onStrand=='+')? 0 : 1;
			purineTrack = duplex[strand];
			if (strand == 1){ // reverse if crick strand
				purineTrack = purineTrack.split("").reverse().join("");
			}
			
			guaninePositions = [];
			numberGuanines = 0;
			for (var gi=0;gi < duplex[strand].length; ++gi){
				if (duplex[strand][gi]=="G"){
					guaninePositions.push(1)
					numberGuanines++;
				} else {
					guaninePositions.push(0)
				}
			}
			
			for (i=0;i<json["aaData"].length;i++){
				json["aaData"][i][0] = Array.repeat(' ', json["aaData"][i][3]).join('') + json["aaData"][i][0]+ Array.repeat(' ', onSize-json["aaData"][i][3]-json["aaData"][i][1]).join('');
		    	// compute overlap
		    	var offStart = json["aaData"][i][3];
		    	var offEnd = offStart+json["aaData"][i][7];
		    	var from = onStart > offStart ? onStart : offStart;
		    	var to = onEnd < offEnd ? onEnd : offEnd;
		    	var overlap = to-from;
				json["aaData"][i][1] = overlap ; 
			}
			// mark the active interval
			targetseq = json["aoColumns"][0]["sTitle"];
			seq = targetseq;
			seq += '<br/>' +Array.repeat(' ', onStart).join('')
			for (i=onStart; i<onEnd;i++){
				if ((onStrand=='+' && (targetseq[i]=='G' || targetseq[i]=='A')) || (onStrand=='-' && (targetseq[i]=='C' || targetseq[i]=='T'))){
					seq += '-';
				} else {
					seq += 'x';
				}
			}
			seq += Array.repeat(' ', onSize - onEnd).join('') ;
			json["aoColumns"][0]["sTitle"] = seq;

			// add details button
			json["fnDrawCallback"]= function ( oSettings ) {
	            $('#offtarget_table tbody tr td.btn_details').each( function (i) {
	            	if ($(this).find("img").is('.btn_detail')) {
	            		return true;
	            	} else {
		            	$(this).replaceWith('<td class="btn_details"><img src="./includes/image/details_open.png" class="btn_detail" /></td>');
	            	}
	            });
	            // also calculate mismatches
	           $.fn.dataTableExt.oApi.fnCalculateMismatches(oSettings);
	        };	 
        
	        // create table
			dTable= $('#offtarget_table').dataTable(json);
			
			// data obtained now start processing
			$('#datanotice').css({opacity: 1.0, visibility: "hidden"}).animate({opacity: 0.0});
			dTable.fnProcessingIndicator( true );
			
			new FixedHeader( dTable ); 
			
			fnCollapseOfftargetCategories(dTable);
			
			// add table footer
			for (i =0;i<$('#offtarget_table thead tr').children('th').length;++i){
				$("#offtarget_table tfoot tr").append('<th style="vertical-align: top; padding-top: 10px"></th>');
			}

   	       	// add filter in the footer
	       	$("#offtarget_table tfoot tr th").each( function ( i ) {
				if (i==1){ // overlap show up to n-3 shorter once by default
		    		this.innerHTML = fnCreateInputSlider( i, dTable.fnGetColumnData(dTable.fnGetTh(i)) );
		    		var extremes = getExtremes(dTable.fnGetColumnData(dTable.fnGetTh(i)) );
		    		fnUpdateInputSlider( dTable, i, dTable.fnGetColumnData(dTable.fnGetTh(i)), parameters["minLength"], extremes[1], 1, Math.max(extremes[0], extremes[1]-3), extremes[1] )	       			
	       		}
		    	else if (i == 2 || i == 4){ // count values
		    		this.innerHTML = fnCreateInputSlider( i, dTable.fnGetColumnData(dTable.fnGetTh(i)) );
		    		var extremes = getExtremes(dTable.fnGetColumnData(dTable.fnGetTh(i)) );
		    		fnUpdateInputSlider( dTable, i, dTable.fnGetColumnData(dTable.fnGetTh(i)), extremes[0], extremes[1], 1, extremes[0], extremes[1] )
		    	} 
		    	else if (i == 5){ // chromatin values
		    		this.innerHTML = fnCreateInputSlider( i, dTable.fnGetColumnData(dTable.fnGetTh(i)) );
		    		var extremes = getExtremes(dTable.fnGetColumnData(dTable.fnGetTh(i)) );
					if (extremes[0] != '-' && extremes[1] != '-'){
		    			fnUpdateInputSlider( dTable, i, dTable.fnGetColumnData(dTable.fnGetTh(i)), extremes[0], extremes[1], (extremes[1]-extremes[0])/25. , extremes[0], extremes[1] );
		    			chromatin = true;
					}
		    	}
		    } );
		    		    
		    // add filtering for invalid rows
			$.fn.dataTableExt.afnFiltering.push(
			    function( oSettings, aData, iDataIndex ) {
			    	var qualifies = true;
			    	// filter with annotation
			    	if (!($.getUrlVar('annotation') == "-" || $.getUrlVar('annotation') == "") && aData[offTargetRegionIds[$.getUrlVar('annotation')]] == 0){
			    		qualifies = false;
			    		return qualifies;
			    	}
			    	
			    	// check against slider values
					for (iColumn=1; iColumn<6; iColumn++){
						if (iColumn == 1 || iColumn == 2 || iColumn == 4 || iColumn == 5 && chromatin ){
							var iDataColumn = dTable.fnGetTh(iColumn);
					    	var iMin = $( "#slider-" + iColumn ).slider( "values", 0 );
					        var iMax = $( "#slider-" + iColumn ).slider( "values", 1 );
					        var iVersion = aData[iDataColumn] == "-" ? 0 : aData[iDataColumn]*1;
					        if ( iMin === "" && iMax === "" ){
					            continue;
					        } else if ( iMin === "" && iVersion < iMax ){
					            continue;
					        } else if ( iMin <= iVersion && "" === iMax ){
					            continue;
					        } else if ( iMin <= iVersion && iVersion <= iMax ){
					            continue;
					        }
					        // violates in at least one column
					        qualifies = false;
					        return qualifies;
						}
			    	}

  			    	return qualifies;
			    }
			);

			// add interaction mode for detail button
			$('#offtarget_table tbody td img').live('click', function () {
			   var nTr = this.parentNode.parentNode;
			   if ( this.src.match('details_close') ){
			       /* This row is already open - close it */
			       this.src = "./includes/image/details_open.png";
			       dTable.fnClose( nTr );
			   } else {
			       /* Open this row */
			       this.src = "./includes/image/details_close.png";
			       dTable.fnOpen( nTr, fnFormatDetails(dTable, nTr), 'details' );
			   }
			} );
			
		    // add tooltip		   
			$('#offtarget_table thead tr th').each( function(i) {
				this.setAttribute( 'title', json["aoColumns"][dTable.fnGetTh(i)]["sToolTip"] );
			} );
			// Apply the tooltips
			$('#offtarget_table thead tr th[title]').tooltip( {
				position:'bottom center',
				slideInSpeed:500,
				opacity:1,
				predelay:2500,
				direction:'up',
				effect:'slide'
			} );
			
			dTable.fnDraw();
			
			dTable.fnProcessingIndicator( false );
			// just to be sure
			$('#datanotice').css('visibility','hidden');

		},
		"error": handleAjaxError // this sets up jQuery to give me errors
		
	}); 
}
