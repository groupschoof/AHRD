package ahrd.model;

import static ahrd.model.ReferenceDescription.tokenizeDescription;
import static ahrd.controller.Settings.getSettings;

import java.util.Iterator;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.search.SearchContentAdapter;

/**
 * Adapter enables memory-friendly parsing of Blast-Result-Files.
 * 
 * see http://www.biojava.org/wiki/BioJava:CookBook:Blast:Echo
 * 
 * On endSubHit (the end of the FIRST HSP) the result is stored in the
 * Query-Protein's BlastResult-Hash.
 * 
 * @author klee, hallab
 */
public class BlastSearchContentAdapter extends SearchContentAdapter {

	private Map<String, Protein> proteinDb;
	private String blastDatabaseName;
	private BlastResult currentBlastResult;
	// Temporary object for collecting the highest scoring HSP:
	private BlastResult bestScoringHSP;
	private String currentProteinAccession;

	public BlastSearchContentAdapter(Map<String, Protein> proteinDb,
			String blastDatabaseName) {
		setProteinDb(proteinDb);
		setBlastDatabaseName(blastDatabaseName);
	}

	public void findHighestScoringHSP() {
		if (this.bestScoringHSP == null
				|| this.currentBlastResult.getBitScore() > this.bestScoringHSP
						.getBitScore()) {
			this.bestScoringHSP = this.currentBlastResult;
			this.currentBlastResult = new BlastResult(this.currentBlastResult
					.getAccession(), this.currentBlastResult.getDescription());
		}
	}

	public static String validateDouble(String input) {
		if (input.startsWith("e") || input.startsWith("E"))
			input = "1" + input;
		return input;
	}

	public boolean passesBlacklist(String blastResultDescriptionLine) {
		boolean passesBlacklist = (blastResultDescriptionLine != null && !blastResultDescriptionLine
				.equals(""));
		for (Iterator<String> i = getSettings().getBlastResultsBlackList(
				getBlastDatabaseName()).iterator(); (i.hasNext() && passesBlacklist);) {
			Pattern p = Pattern.compile(i.next());
			Matcher m = p.matcher(blastResultDescriptionLine);
			passesBlacklist = !m.find();
		}
		return passesBlacklist;
	}

	public String filter(String blastResultDescriptionLine) {
		String filteredDescLine = blastResultDescriptionLine;
		for (Iterator<String> i = getSettings().getBlastResultsFilter(
				getBlastDatabaseName()).iterator(); i.hasNext();) {
			Pattern p = Pattern.compile(i.next());
			// Replace with whitespace, so word-boundaries are kept up
			filteredDescLine = p.matcher(filteredDescLine).replaceAll(" ");
		}
		// Condense multiple whitespaces into one and trim the description-line:
		filteredDescLine = filteredDescLine.replaceAll("\\s{2,}", " ").trim();
		return filteredDescLine;
	}

	public void startHit() {
		// initialize Hits
		this.currentBlastResult = new BlastResult(getBlastDatabaseName());
		this.bestScoringHSP = null;
	}

	public void endSubHit() {
		findHighestScoringHSP();
	}

	/**
	 * Adds the currently best scoring BlastHit (with alignment-values from it's
	 * best scoring HSP) to the Protein's BlastResults, if and only if this
	 * BlastResult passes blacklist and filtering.
	 * 
	 * @Note: If requested, remembers the unchanged BlastHit for AHRD-training
	 *        and evaluation - in this case the BlastHit is not filtered and not
	 *        passed through the blacklist.
	 */
	public void endHit() {
		if (this.bestScoringHSP != null) {
			if (getProteinDb().containsKey(currentProteinAccession)) {
				Protein protein = getProteinDb().get(currentProteinAccession);
				// For Training-Purposes:
				if (getSettings().getWriteBestBlastHitsToOutput()) {
					// Of course we do have to treat this best-blast-hit
					// differently than the further to process one below, so
					// clone:
					BlastResult theClone = this.bestScoringHSP.clone();
					// Pass best Blast-Hit's Description through filter:
					theClone.setDescription(filter(theClone.getDescription()));
					// Tokenize without filtering tokens through the Blacklist:
					theClone.setTokens(tokenizeDescription(theClone
							.getDescription()));
					protein.getEvaluationScoreCalculator()
							.addUnchangedBlastResult(getBlastDatabaseName(),
									theClone);
				}
				if (passesBlacklist(this.bestScoringHSP.getDescription())) {
					// Pass bestScoringHSP through filter:
					this.bestScoringHSP
							.setDescription(filter(this.bestScoringHSP
									.getDescription()));
					// Tokenize the filtered Description-Line:
					this.bestScoringHSP.tokenize();
					// Pass bestScoringHSP through Blacklist and add it, if it
					// is still valid:
					if (this.bestScoringHSP.isValid()) {
						// Generate the pattern of the Description's unique
						// tokens:
						this.bestScoringHSP.patternize();
						// Adds the BlastResult to the Protein's set and
						// measures the cumulative and total scores later needed
						// to calculate the Token-Scores:
						protein.addBlastResult(this.bestScoringHSP);
					}
				}
			} else {
				System.err
						.println("ERROR: Could not find Protein for Accession '"
								+ currentProteinAccession + "'.");
			}
		}
	}

	/**
	 * Handles Blast-Hits
	 */
	public void addHitProperty(Object key, Object val) {
		String theKey = key.toString();
		String theValue = val.toString().trim();
		if (theKey.equals("subjectId")) {
			currentBlastResult.setAccession(theValue);
		} else if (theKey.equals("subjectDescription")) {
			currentBlastResult.setDescription(theValue);
		}
	}

	/**
	 * Handles Blast-HSPs
	 */
	public void addSubHitProperty(Object key, Object val) {
		String theKey = key.toString();
		String theValue = val.toString().trim();
		if (theKey.equals("score")) {
			currentBlastResult.setBitScore(Double.parseDouble(theValue));
		} else if (theKey.equals("expectValue")) {
			currentBlastResult.setEValue(Double
					.valueOf(validateDouble(theValue)));
		} else if (theKey.equals("subjectSequenceStart")) {
			currentBlastResult.setStart(Integer.parseInt(theValue));
		} else if (theKey.equals("subjectSequenceEnd")) {
			currentBlastResult.setEnd(Integer.parseInt(theValue));
		}
	}

	public void setQueryID(String queryID) {
		currentProteinAccession = queryID.trim();
	}

	public Map<String, Protein> getProteinDb() {
		return proteinDb;
	}

	public void setProteinDb(Map<String, Protein> proteinDb) {
		this.proteinDb = proteinDb;
	}

	public String getBlastDatabaseName() {
		return blastDatabaseName;
	}

	public void setBlastDatabaseName(String blastDatabaseName) {
		this.blastDatabaseName = blastDatabaseName;
	}

}
