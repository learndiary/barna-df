package fbi.commons.tools;

import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;

import fbi.commons.flux.FluxTool;
import fbi.commons.options.HelpPrinter;

@Cli(name = "subsetter", description = "Extract a random subset of lines from a file")
public class Subsetter implements FluxTool<Void> {

	
	
	
	@Override
	public Void call() throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean validateParameters(HelpPrinter printer,
			ArgumentProcessor toolArguments) {
		// Log.error
		return true;
	}
	
}
