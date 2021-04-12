package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import ahrd.controller.AHRD;
import ahrd.model.BlastResult;
import ahrd.model.GOterm;
import ahrd.model.Protein;

public class TsvOutputWriter extends OutputWriter {

	protected BufferedWriter hrdScoresWriter;

	public TsvOutputWriter(Collection<Protein> proteins) {
		super(proteins);
	}

	public void writeOutput() throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(getSettings().getPathToOutput()));
		if (getSettings().doWriteHRDScoresToOutput())
			writeHRDScoresOutputHeader();

		// Header:
		bw.write("# AHRD-Version " + AHRD.VERSION + "\n");
		bw.write("\n");
		// Column-Names:
		bw.write(ahrdColumnNames());
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGoSlimFile()) {
			bw.write("\tGO-Slim-Annotation");
		}
		bw.write("\n");

		for (Protein prot : getProteins()) {
			// Generate the Human Readable Description:
			String csvRow = buildDescriptionLine(prot, "\t");
			if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGoSlimFile()) {
				csvRow += "\t" + combineGoTermsToString(prot.getGoSlimTerms());
			}
			
			// Write row to CSV:
			csvRow += "\n";
			bw.write(csvRow);

			// If AHRD is requested to write out the AHRD-Score of each
			// BlastHit's Description, do so into another file:
			if (getSettings().doWriteHRDScoresToOutput())
				writeHrdScoresOutput(prot);
		}

		// CLEAN UP:
		bw.close();
		if (getSettings().doWriteHRDScoresToOutput())
			this.hrdScoresWriter.close();
	}
	
	/** 
	 * @return String - Tab separated column names for the basic AHRD-output
	 */
	public String ahrdColumnNames() {
		return "Protein-Accession\tBlast-Hit-Accession\tAHRD-Quality-Code\tHuman-Readable-Description\tGene-Ontology-Term";
	}

	/**
	 * AHRD can be requested to log all final AHRD-Scores of each BlastHit's
	 * Description. This is needed for fitting a generalized extreme value
	 * distribution to the AHRD-Scores. The fitted gevd can later be used to
	 * calculate P-Values for the AHRD-Scores assigned to each BlastHit's
	 * Description. This method initializes the TsvOutputWriter for the above
	 * scores.
	 */
	public void writeHRDScoresOutputHeader() throws IOException {
		// Initialize TsvOutputWriter:
		hrdScoresWriter = new BufferedWriter(new FileWriter(getSettings().getPathToHRDScoresOutput()));
		hrdScoresWriter.write("Protein-Accesion\tBlast-Hit-Accession\tAHRD-Score\n");
	}

	/**
	 * AHRD can be requested to log all final AHRD-Scores of each BlastHit's
	 * Description. This is needed for fitting a generalized extreme value
	 * distribution to the AHRD-Scores. The fitted gevd can later be used to
	 * calculate P-Values for the AHRD-Scores assigned to each BlastHit's
	 * Description. This method writes the above scores for the argument
	 * protein.
	 * 
	 * @throws IOException
	 */
	public void writeHrdScoresOutput(Protein prot) throws IOException {
		for (String blastDatabaseName : prot.getBlastResults().keySet()) {
			for (BlastResult br : prot.getBlastResults().get(blastDatabaseName)) {
				this.hrdScoresWriter
						.write(prot.getAccession() + "\t" + br.getAccession() + "\t" + br.getDescriptionScore() + "\n");
			}
		}
	}
	
	public String combineGoTermsToString(Set<GOterm> gos) {
		return combineGoTermsToString(gos, ", ");
	}

	public String combineGoTermsToString(Set<GOterm> gos, String separator) {
		String goLine = "";
		if (gos != null) {
			List<String> sortedGos = new ArrayList<String>();
			for (Iterator<GOterm> goTermIter = gos.iterator(); goTermIter.hasNext();) {
				sortedGos.add(goTermIter.next().getAccession());
			}
			Collections.sort(sortedGos);
			for (Iterator<String> iter = sortedGos.iterator(); iter.hasNext();) {
				goLine += iter.next();
				if (iter.hasNext())
					goLine += separator;
			}

		}
		return goLine;
	}
	
}
