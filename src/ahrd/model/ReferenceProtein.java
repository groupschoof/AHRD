package ahrd.model;

import static ahrd.controller.Settings.FASTA_PROTEIN_HEADER_ACCESSION_GROUP_NAME;
import static ahrd.controller.Settings.FASTA_PROTEIN_HEADER_DESCRIPTION_GROUP_NAME;
import static ahrd.controller.Settings.GO_TERM_GROUP_NAME;
import static ahrd.controller.Settings.SHORT_ACCESSION_GROUP_NAME;
import static ahrd.controller.Settings.getSettings;
import static ahrd.model.AhrdDb.getReferenceProteinDAO;
import static com.sleepycat.persist.model.Relationship.ONE_TO_ONE;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.sleepycat.persist.model.Entity;
import com.sleepycat.persist.model.PrimaryKey;
import com.sleepycat.persist.model.SecondaryKey;

import ahrd.controller.Settings;

@Entity
public class ReferenceProtein {

	/**
	 * Store (persist) the ReferenceProtein.
	 * 
	 * @param rp
	 */
	public static void put(ReferenceProtein rp) {
		getReferenceProteinDAO().byAccession.put(rp);
	}

	/**
	 * Query the persistence layer for the ReferenceProtein belonging to the
	 * argument accession.
	 * 
	 * @param referenceProteinAccession
	 *            can be either the short or the normal accession
	 * @return ReferenceProtein matching the accession
	 */
	public static ReferenceProtein get(String referenceProteinAccession) {
		ReferenceProtein rp = getReferenceProteinDAO().byAccession.get(referenceProteinAccession);
		if (rp == null) {
			rp = getReferenceProteinDAO().byShortAccession.get(referenceProteinAccession);
		}
		return rp;
	}

	/**
	 * Extracts from the provided protein database in FASTA format the length
	 * and human readable descriptions of those Proteins that are Hits
	 * (Subjects) in the argument blastResults. Each time such a Hit is found
	 * the mentioned measurements are set in the respective instance of
	 * BlastResult and subsequently the method 'Protein.addBlastResult' is
	 * invoked.
	 * 
	 * @param blastDbName
	 * @throws IOException
	 */
	public static void parseBlastDatabase(String blastDbName) throws IOException {
		// Parse line by line FASTA Blast search DB. Extract Subject Lengths and
		// Subject HRDs.
		BufferedReader fastaIn = null;
		try {
			fastaIn = new BufferedReader(new FileReader(getSettings().getPathToBlastDatabase(blastDbName)));
			String str, hrd = new String();
			String acc = "";
			Long hitAALength = new Long(0);
			boolean hit = false;
			while ((str = fastaIn.readLine()) != null) {
				if (str.startsWith(">")) {
					// Finished reading in the original Fasta-Entry of a
					// Blast-Hit? If so, process it:
					if (hit) {
						put(new ReferenceProtein(acc, hrd, hitAALength, blastDbName));
						// Clean up to enable processing the next Hit
						hitAALength = new Long(0);
						// Note, that the boolean 'hit' will be set in the
						// following If-Else-Block.
					}
					// Process the current Fasta-Header-Line:
					Matcher m = getSettings().getFastaHeaderRegex(blastDbName).matcher(str);
					if (!m.matches()) {
						// Provided REGEX to parse FASTA header does not work in
						// this case:
						System.err.println("WARNING: FASTA header line\n" + str.trim()
								+ "\ndoes not match provided regular expression\n"
								+ getSettings().getFastaHeaderRegex(blastDbName).toString()
								+ "\n. The header and the following entry, including possibly respective matching BLAST Hits, are ignored and discarded.\n"
								+ "To fix this, please use - Blast database specific - parameter "
								+ Settings.FASTA_HEADER_REGEX_KEY
								+ " to provide a regular expression that matches ALL FASTA headers in Blast database '"
								+ blastDbName + "'.");
					} else {
						// Found the next Reference Protein:
						acc = m.group(FASTA_PROTEIN_HEADER_ACCESSION_GROUP_NAME).trim();
						hrd = m.group(FASTA_PROTEIN_HEADER_DESCRIPTION_GROUP_NAME).trim();
						// Following lines, until the next header, contain
						// information to be collected:
						hit = true;
					}
				} else if (hit) {
					// Process non header-line, if and only if, we are reading
					// the sequence of a Blast-Hit:
					hitAALength += str.trim().length();
				}
			}
			// Was the last read FASTA entry a Blast-Hit? If so, it needs
			// processing:
			if (hit)
				put(new ReferenceProtein(acc, hrd, hitAALength, blastDbName));
		} finally {
			fastaIn.close();
		}
	}

	/**
	 * Parses the tabular reference Gene Ontology term annotations (GOA) for
	 * proteins in the searched Blast databases. These GOAs will then be used to
	 * annotate the query proteins with GO terms. The important restriction is,
	 * that only those reference GO annotations will be extracted that match one
	 * of the BlastResults found in the respective Blast searches
	 * <code>uniqueShortAccessions</code>.
	 * 
	 * @param Set
	 *            <String> uniqueShortAccessions - The unique short accessions
	 *            of BlastResults found in the respective Blast searches
	 * @throws IOException
	 */
	public static void parseReferenceGoAnnotations() throws IOException {
		BufferedReader goaIn = null;
		try {
			goaIn = new BufferedReader(new FileReader(getSettings().getPathToGeneOntologyResults()));
			Pattern p = getSettings().getReferenceGoRegex();
			String line, shortAcc, goTerm = "";
			while ((line = goaIn.readLine()) != null) {
				Matcher m = p.matcher(line);
				if (m.find()) {
					shortAcc = m.group(SHORT_ACCESSION_GROUP_NAME);
					goTerm = m.group(GO_TERM_GROUP_NAME);
					addGoAnnotation(shortAcc, goTerm);
				}
			}
		} finally {
			goaIn.close();
		}
	}

	/**
	 * Adds the Gene Ontology term <code>goTerm</code> to the Set of the
	 * BlastResult's GO term annotations. In this, the BlastResult is identified
	 * by its short accession <code>brShortAccession</code>.
	 * 
	 * @param brShortAccession
	 * @param goTerm
	 */
	protected static void addGoAnnotation(String brShortAccession, String goTerm) {
		ReferenceProtein rp = getReferenceProteinDAO().byShortAccession.get(brShortAccession);
		if (null != rp) {
			// Add Go Term
			rp.getGoTerms().add(goTerm);
			// Update Database-Entry:
			getReferenceProteinDAO().byAccession.put(rp);
		} else
			System.err.println("WARNING: Could not find Reference-Protein for SHORT-ACCESSION '" + brShortAccession
					+ "' in the AHRD Database. But found GO-Terms for it. Discarding them...");
	}

	@PrimaryKey
	private String accession;
	@SecondaryKey(relate = ONE_TO_ONE)
	private String shortAccession;
	private String hrd;
	private Long sequenceLength;
	String blastDatabaseName;
	private Set<String> goTerms = new HashSet<String>();

	/**
	 * For deserialization
	 */
	public ReferenceProtein() {
	}

	/**
	 * Constructor setting the fields
	 * 
	 * @param accession
	 * @param hrd
	 * @param sequenceLength
	 */
	public ReferenceProtein(String accession, String hrd, Long sequenceLength, String blastDatabaseName) {
		super();
		setAccession(accession);
		setHrd(hrd);
		setSequenceLength(sequenceLength);
		setBlastDatabaseName(blastDatabaseName);
		// Extract the short accession from the long one:
		initShortAccession();
	}

	/**
	 * Extracts from the possibly longer Accession the shorter one which is used
	 * in the reference Gene Ontology Annotation (GOA) file. This is required
	 * due to the fact that UniprotKB uses short accessions in the GOA files but
	 * provides long accessions in the Blast Databases. A very unfortunate lack
	 * of standardization, indeed!
	 */
	public void initShortAccession() {
		Pattern p = getSettings().getShortAccessionRegex(getBlastDatabaseName());
		Matcher m = p.matcher(getAccession());
		setShortAccession(getAccession());
		if (!m.find()) {
			System.err.println("WARNING: Regular Expression '" + p.toString()
					+ "' does NOT match - using pattern.find(...) - Blast Hit Accession '" + getAccession()
					+ "' - continuing with the original accession. This might lead to unrecognized reference GO annotations!");
		} else {
			setShortAccession(m.group(SHORT_ACCESSION_GROUP_NAME));
		}
	}

	public String getShortAccession() {
		return this.shortAccession;
	}

	public void setShortAccession(String shortAccession) {
		this.shortAccession = shortAccession;
	}

	public String getAccession() {
		return accession;
	}

	public void setAccession(String accession) {
		this.accession = accession;
	}

	public String getHrd() {
		return hrd;
	}

	public void setHrd(String hrd) {
		this.hrd = hrd;
	}

	public Long getSequenceLength() {
		return sequenceLength;
	}

	public void setSequenceLength(Long sequenceLength) {
		this.sequenceLength = sequenceLength;
	}

	public Set<String> getGoTerms() {
		return goTerms;
	}

	public void setGoTerms(Set<String> goTerms) {
		this.goTerms = goTerms;
	}

	public String getBlastDatabaseName() {
		return blastDatabaseName;
	}

	public void setBlastDatabaseName(String blastDatabaseName) {
		this.blastDatabaseName = blastDatabaseName;
	}
}
