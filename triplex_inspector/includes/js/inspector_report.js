var rTable = null;
var oTable = null;
var thisUrlVars; 

// open new window and show detail for given region
function showDetail(rId){
	window.open('inspector_detail.html?rId='+rId+'&region=&annotation=', '_blank');
}

function addFilter(){
	/* preselect annotation  */
	$('#filter_table tfoot tr th:eq(8) select option').filter(function() {
		//may want to use $.trim in here
		return $(this).text().trim() == "-"; 
	}).attr('selected', true).change();

	/* Add events */
    $('#onregion_table tbody tr').live('click', function () {
        var nTds = $('td', this);
        var regionId = $(nTds[0]).text();
        $('#filter_table tfoot tr th:eq(0) option:selected').attr('selected', false);
        $('#filter_table tfoot tr th:eq(0) select option').filter(function() {
		    //may want to use $.trim in here
		    return $(this).text() == regionId; 
		}).attr('selected', true).change();
    } );
    
	$("#onregion_table tbody").click(function(event) {
        $(rTable.fnSettings().aoData).each(function (){
            $(this.nTr).removeClass('row_selected');
        });
        $(event.target.parentNode).addClass('row_selected');
    });
    
    $('#filter_table tbody tr').live('click', function () {
        var nTds = $('td', this);
        var regionId = $(nTds[0]).text();
        var region = $(nTds[1]).text();
        var annotation = $(nTds[8]).text();
		window.open('inspector_detail.html?rId='+regionId+'&region='+region+'&annotation='+annotation, '_blank');
    } );
    
	$("#filter_table tbody").click(function(event) {
        $(oTable.fnSettings().aoData).each(function (){
            $(this.nTr).removeClass('row_selected');
        });
        $(event.target.parentNode).addClass('row_selected');
    });	
}
	
function loadPrimaryTargets(){
	// load all primary targets
	$.ajax({
		url: 'json/'+submatches_file,
		async: false,
		dataType: 'json',
		success: function (json) {
//			json["sDom"] = 'C<"clear">lfrtip'; // deactivate due to bug for the time
			json["aLengthMenu"] = [[10, 25, 50, -1], [10, 25, 50, "All"]];
//			json["bStateSave"] = true; // TODO conflicts with filters at table bottom
			json["fnRowCallback"] = function( nRow, aData, iDisplayIndex ) {
				// indicate shorter off-target categories that have not been processed due to parameter setting
				if (aData[10] == 0 && aData[2]-1 < parameters["minLength"]){
					jQuery('td:eq(10)', nRow).html( '-' );
				}
				if (aData[11] == 0 && aData[2]-2 < parameters["minLength"]){
					jQuery('td:eq(11)', nRow).html( '-' );
				}
				return nRow;
			};
			
			oTable= $('#filter_table').dataTable(json);
			$('#datanotice').css({opacity: 1.0, visibility: "hidden"}).animate({opacity: 0.0});
			oTable.fnProcessingIndicator( true );
			
			new FixedHeader( oTable ); 
			
	       	// add filter in the footer
	       	$("#filter_table tfoot th").each( function ( i ) {
		    	if (i == 0 || i == 8){ // categorical values
			        this.innerHTML = fnCreateSelect( oTable.fnGetColumnData(i) );
			        $('select', this).change( function () {
			            oTable.fnFilter( $(this).val(), i, false, false );
			        } );
		    	} else if (i == 2 || i == 4 || i == 5 || i == 6 || i == 9){ // count values
		    		this.innerHTML = fnCreateInputSlider( i, oTable.fnGetColumnData(i) );
		    		var extremes = getExtremes(oTable.fnGetColumnData(i) );
		    		fnUpdateInputSlider( oTable, i, oTable.fnGetColumnData(i), extremes[0], extremes[1], 1, extremes[0], extremes[1] );
		    	} else if (i == 3 || i == 7 ){ // ratio values
		    		this.innerHTML = fnCreateInputSlider( i, oTable.fnGetColumnData(i));
		    		var extremes = getExtremes(oTable.fnGetColumnData(i) );
		    		fnUpdateInputSlider( oTable, i, oTable.fnGetColumnData(i), extremes[0], extremes[1], 0.01, extremes[0], extremes[1] );
		    	} 
		    } );
		    
		    $('#filter_table').dataTable().fnSettings().aoDrawCallback.push( {
			    "fn": function () {
				    	$("#filter_table tbody tr").each( function (i){
							var aData = oTable.fnGetData( this );
					
							if (aData[9] <= 10){
								$(this).addClass("zeroOffs");
							} else if (aData[9] <= 100) {
								$(this).addClass("tensOffs");
							} else if (aData[9] <= 1000) {
								$(this).addClass("hundredsOffs");
							} else {
								$(this).addClass("thousandsOffs");
							}
						});
			    	}
			} );
		    
		    
		    // add custom filtering for columns 3-7
			$.fn.dataTableExt.afnFiltering.push(
			    function( oSettings, aData, iDataIndex ) {
			    	var qualifies = true;
			    	for (var iColumn=2; iColumn<=9; iColumn++){
			    		if (iColumn==8) continue; // skip annotation column
				    	var iMin = $( "#slider-" + iColumn ).slider( "values", 0 );
				        var iMax = $( "#slider-" + iColumn ).slider( "values", 1 );
				        
				        var iVersion = aData[iColumn] == "-" ? 0 : aData[iColumn]*1;
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
			    	}
			    	return qualifies;
			    }
			);
			
		    // add tooltip		   
			$('#filter_table thead tr:last-child th').each( function(i) {
				this.setAttribute( 'title', json["aoColumns"][oTable.fnGetTh(i)]["sToolTip"] );
			} );
			// Apply the tooltips
			$('#filter_table thead tr:last-child th[title]').tooltip( {
				position:'bottom center',
				slideInSpeed:500,
				opacity:1,
				predelay:2000,
				direction:'up',
				effect:'slide'
			} );
			
			$('#filter_table').dataTable().fnDraw();
			oTable.fnProcessingIndicator( false );

			// just to be sure
			$('#datanotice').css('visibility','hidden');  
			
			addFilter();
		}
	});	
}	
		
$(document).ready(function(){
	// load primary target regions
	$.ajax({
		url: 'json/primary_target_regions.json',
		async: false,
		dataType: 'json',
		success: function (json) {
			// json works, hide warning
			$('#browsererror').css('visibility','hidden');
			$('#datanotice').css({opacity: 0.0, visibility: "visible"}).animate({opacity: 1.0});

			rTable= $('#onregion_table').dataTable(json);
			
		    // add tooltip		   
			$('#onregion_table thead tr th').each( function(i) {
				this.setAttribute( 'title', json["aoColumns"][rTable.fnGetTh(i)]["sToolTip"] );
			} );
			// Apply the tooltips
			$('#onregion_table thead tr th[title]').tooltip( {
				position:'bottom center',
				slideInSpeed:500,
				opacity:1,
				predelay:2000,
				direction:'up',
				effect:'slide'
			} );		
			$('#onregion_table').dataTable().fnDraw();
		}
	});

	setTimeout("loadPrimaryTargets()",250);
  
});
