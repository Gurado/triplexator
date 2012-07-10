$(window).load(function() {
	var width2height = $(window).width() / $(window).height();
	var colSpots = Math.ceil(Math.sqrt($('img').length));
	var scale = Math.min($(window).height(), $(window).width()/colSpots) / $('img').first().width();
	var scaled_imgage_width = scale*$('img').first().width();
	var row = 0, col = 0;
	for (var i=0; i<$('img').length; i++){
		img = $('img:eq('+i+')');
		// scale
		img.css('-webkit-transform', 'scale('+scale+','+scale+')');
		img.css('-moz-transform', 'scale('+scale+','+scale+')');
		img.css('-ms-transform', 'rscale('+scale+','+scale+')');
		img.css('transform', 'scale('+scale+','+scale+')');
		// position
		img.css('position', 'absolute');
		img.css('top', scaled_imgage_width * row);
		img.css('left', scaled_imgage_width * col);
		col++;
		if (col >= colSpots){
			col = 0;
			row++;
		}
	}
});
