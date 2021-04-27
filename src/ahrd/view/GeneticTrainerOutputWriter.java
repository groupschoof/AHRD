package ahrd.view;

import java.io.IOException;

import ahrd.controller.Parameters;

public class GeneticTrainerOutputWriter extends TrainerOutputWriter {

	public GeneticTrainerOutputWriter(String logPath) throws IOException {
		super(logPath);
	}

	@Override
	public String generateHeader(boolean isFinalOutput, Parameters p) {
		String hdr = "Generation" + separator;
		if (isFinalOutput)
			hdr += "Average Maximum-Evaluation-Score" + separator;
		hdr += "Average Evaluation-Score(F-Score)" + separator;
		if (!isFinalOutput)
			hdr += "Diff-to-last-Generation" + separator + "Origin" + separator;
		hdr += "Average Precision" + separator;
		hdr += "Average Recall" + separator;
		hdr += p.buildHeaderForOutput(separator);
		hdr += "\n";
		return hdr;
	}
	
}
