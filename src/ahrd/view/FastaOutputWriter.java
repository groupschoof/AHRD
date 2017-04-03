package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import ahrd.model.Protein;

public class FastaOutputWriter extends OutputWriter {

	public FastaOutputWriter(Collection<Protein> proteins) {
		super(proteins);
	}

	public void writeOutput() throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(getSettings()
				.getPathToOutput()));

		for (Protein prot : getProteins()) {
			// Write Fasta-Header
			bw.write(">" + buildDescriptionLine(prot, " ") + "\n");
			// Append AA-Sequence
			bw.write(prot.getSequence() + "\n");
		}

		bw.close();
	}
}
