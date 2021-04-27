package ahrd.view;

import java.io.IOException;

import ahrd.controller.Parameters;

public class SimulatedAnnealingTrainerOutputWriter extends TrainerOutputWriter {

	public SimulatedAnnealingTrainerOutputWriter(String logPath) throws IOException {
		super(logPath);
	}
	
	@Override
	public String generateHeader(boolean isFinalOutput, Parameters p) {
		String hdr = "Temperature" + separator;
		if (isFinalOutput)
			hdr += "Average Maximum-Evaluation-Score" + separator;
		hdr += "Average Evaluation-Score(F-Score)" + separator;
		if (!isFinalOutput)
			hdr += "Diff-to-curr-Accepted" + separator + "Accepted"  + separator;
		hdr += "Average Precision" + separator;
		hdr += "Average Recall" + separator;
		hdr += p.buildHeaderForOutput(separator);
		hdr += "\n";
		return hdr;
	}

}
