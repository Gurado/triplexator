var thisUrlVars; 

// get the parameters passed to this url
$.extend({
	getUrlVars: function(){
		var vars = [], hash;
		var hashes = window.location.href.slice(window.location.href.indexOf('?') + 1).split('&');
		for(var i = 0; i < hashes.length; i++){
			hash = hashes[i].split('=');
			vars.push(hash[0]);
			vars[hash[0]] = hash[1];
		}
		return vars;
	},
	getUrlVar: function(name){
		return $.getUrlVars()[name];
	}
});

// enables to generate spacers
Array.repeat= function(val, len){
    for(var i= len, a= []; i--; ) a[i]= val;
    return a;
}

// strip whitespaces from string
// " id".trim() === "id"
if(typeof(String.prototype.trim) === "undefined")
{
    String.prototype.trim = function() 
    {
        return String(this).replace(/^\s+|\s+$/g, '');
    };
}

// shitch processing notice
jQuery.fn.dataTableExt.oApi.fnProcessingIndicator = function ( oSettings, onoff ) {
    if ( typeof( onoff ) == 'undefined' ) {
        onoff = true;
    }
    this.oApi._fnProcessingDisplay( oSettings, onoff );
};

// maps the table header of visible columns back to the aaData/aoColumns index
$.fn.dataTableExt.oApi.fnGetTh  = function ( oSettings, iTh ){
	var counter = -1
	var visIndex = 0
	for (visIndex=0;visIndex<oSettings.aoColumns.length;++visIndex){
		if (oSettings.aoColumns[visIndex]["bVisible"]){
			++counter;
		} 
		if (counter == iTh){
			break;
		} 
	}
	return visIndex;
}

// create input slider
function fnCreateInputSlider( label, aData ){
	var r = '<div id="slider-' + label + '" class="ui-slider ui-slider-horizontal ui-widget ui-widget-content ui-corner-all"></div>', i, iLen=aData.length;
	r += '<input type="text" value="range" id="slider-info-' + label + '" style="text-align:center; border:0; color:#000000; background-color: transparent;" />'
	return r;
}

// return the minimum and the maximum value of the data
function getExtremes(aData){
	var min_v = null;
	var max_v = null;
	iLen=aData.length;
	if (iLen > 0){
		min_v = aData[0]
		max_v = aData[0]
		for (i = 0; i < iLen; i++) {
				if (min_v == null || +aData[i] < min_v){
					min_v = +aData[i]	
				}
				if (max_v == null || +aData[i] > max_v){
					max_v = +aData[i];
				}
	    }
	}
	return [min_v, max_v];
}

// on change update a slider values and filter table 
function fnUpdateInputSlider( oTable, label, aData, min_v, max_v, step_size, value1, value2 ){
	$( "#slider-" + label).slider({
			range: true,
			min: eval(min_v),
			max: eval(max_v),
			step: step_size,
			values: [ value1, value2 ],
			slide: function( event, ui ) {
				$( "#slider-info-"+label ).val( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
			},
			stop: function() { oTable.fnDraw(); }
		});
	$( "#slider-info-"+label ).val( $( "#slider-" + label ).slider( "values", 0 ) +
		" - " + $( "#slider-" + label ).slider( "values", 1 ) );
}

/* Get the rows which are currently selected */
function fnGetSelected( oTableLocal ){
    var aReturn = new Array();
    var aTrs = oTableLocal.fnGetNodes();     
    for ( var i=0 ; i<aTrs.length ; i++ ){
        if ( $(aTrs[i]).hasClass('row_selected') )
        {
            aReturn.push( aTrs[i] );
        }
    }
    return aReturn;
}

