$(document).ready(function() {
  
  Shiny.addCustomMessageHandler("version_checkHandler",
    function(version) {
      var version_req = "3.1.3";
      var t= version.localeCompare(version_req);
      if (t == -1) window.alert("You are running R version " + version + ". R version 3.2.2 or more recent is required. The software may not run properly with your version.");
    }
  );

	$("#file").click(function() {
		window.alert("Please make sure that you save necessary summaries before you change data.");
	});
  
	$("#file").change(function() {
		Shiny.onInputChange("error_bars", false);
	});

});
