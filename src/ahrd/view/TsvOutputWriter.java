package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ahrd.controller.AHRD;
import ahrd.controller.Settings;
import ahrd.model.BlastResult;
import ahrd.model.GOterm;
import ahrd.model.Protein;

public class TsvOutputWriter extends OutputWriter {

	protected BufferedWriter hrdScoresWriter;
	protected List<String> goCentricTerms = new LinkedList<String>();

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
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGoTermCentricTermsFile()) {
			bw.write(this.goTermCentricColumnNames());
		}
		bw.write("\n");

		for (Protein prot : getProteins()) {
			// Generate the Human Readable Description:
			String csvRow = buildDescriptionLine(prot, "\t");
			// If requested write GoSlimTerms
			if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGoSlimFile()) {
				csvRow += "\t" + combineGoTermsToString(prot.getGoSlimTerms());
			}
			// If requested write GOterm centric protein-term association confidences 
			if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGoTermCentricTermsFile()) {
				for (String termAcc : goCentricTerms) {
					csvRow += "\t" + FRMT.format(prot.getGoCentricTermConfidences().get(termAcc));				
				}
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
		return "Protein-Accession\tBlast-Hit-Accession\tAHRD-Quality-Code\tHuman-Readable-Description\tInterpro-ID (Description)\tGene-Ontology-Term";
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
	
	/**
	 * AHRD can be requested to annotate proteins with go terms in a term centric way.
	 * This behavior can be triggered by suppling the go term centric terms in a separate file, 
	 * in addition to suppling all the information needed to annotate proteins with go terms in the traditional way.
	 * For every go term centric term a column will be put in AHRDs output with association confidence scores with each protein.
	 * This method creates the names for these columns
	 * @return A string with the column names for the go term centric terms
	 * @throws IOException
	 */
	protected String goTermCentricColumnNames() throws IOException {
		// Load set of GO centric terms
		for (String goCentricTermFileEntry : getSettings().getGoTermCentricTerms()) {
			Pattern p = Settings.GO_TERM_CENTRIC_TERMS_FILE_GOTERM_REGEX;
			Matcher m = p.matcher(goCentricTermFileEntry);
			if (m.find()) {
				String termAcc = m.group("goTerm");
				goCentricTerms.add(termAcc);
			}
		}
		String colNames = null; 
		for (String termAcc : goCentricTerms) {
			colNames = colNames + "\t" + termAcc;				
		}
		return colNames;
	}
	
	public String combineGoTermsToString(Set<GOterm> gos) {
		return combineGoTermsToString(gos, ", ");
	}

	public String combineGoTermsToString(Set<GOterm> gos, String seperator) {
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
					goLine += seperator;
			}

		}
		return goLine;
	}
	
}
