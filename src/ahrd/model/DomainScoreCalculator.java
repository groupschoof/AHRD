package ahrd.model;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * We estimate protein similarity based on the formulae given in the article
 * "Protein comparison at the domain architecture level" by Lee and Lee
 * (http://www.biomedcentral.com/1471-2105/10/S15/S5). The more similar a
 * described protein is to the query the more appropriate its description. Hence
 * we introduce a new Domain-Score on the level of BlastResults to reflect this
 * concept.
 * 
 * @author hallab, bangalore, klee
 */
public class DomainScoreCalculator {

	/**
	 * Result from parsing SIMAP's concatonated and preprocessed feature-files.
	 * Preprocessing involves substitution of SIMAP-Hashes with the original
	 * Protein-Accessions.
	 * 
	 * @note: See awk-scripts in directory helper_scripts.
	 */
	private static Map<String, Set<String>> blastResultAccessionsToInterproIds;

	public static Map<String, Set<String>> getBlastResultAccessionsToInterproIds() {
		return blastResultAccessionsToInterproIds;
	}

	public static void setBlastResultAccessionsToInterproIds(
			Map<String, Set<String>> blastResultAccessionsToInterproIds) {
		DomainScoreCalculator.blastResultAccessionsToInterproIds = blastResultAccessionsToInterproIds;
	}

	public static void initializeBlastResultAccessionsToInterproIds()
			throws IOException {
		blastResultAccessionsToInterproIds = new HashMap<String, Set<String>>();
		BufferedReader reader = new BufferedReader(new InputStreamReader(
				new FileInputStream(getSettings()
						.getPathToInterproResults4BlastHits())));
		String accession = "";
		String interproId = "";
		for (String line; (line = reader.readLine()) != null;) {
			// Parse lines exactly like in method
			// InterproResult.parseInterproResult(Map<String, Protein>)
			// You should extract two values: Protein-Accession and Interpro-ID
			// put those results in Map 'blastResultAccessionsToInterproIds'
			if (!getBlastResultAccessionsToInterproIds().containsKey(accession)) {
				getBlastResultAccessionsToInterproIds().put(accession,
						new HashSet<String>());
			}
			getBlastResultAccessionsToInterproIds().get(accession).add(
					interproId);
		}
	}

	private Protein protein;

	public static Double calculateDomainScore(Protein prot) {
		// 1.) Construct the template for the Vector of domain-weights!
		// Copy the Protein's domains into a new Set:
		Set<InterproResult> domains = new HashSet<InterproResult>(
				prot.getInterproResults());
		// Add all BlastResults' domains to above Set, but discard double
		// entries:
		for (String blastDb : prot.getBlastResults().keySet()) {
			for (BlastResult br : prot.getBlastResults().get(blastDb)) {

			}
		}
		return 0.0;
	}

	/**
	 * Gathers all distinct InterproIDs assigned to the Protein and its
	 * BlastResults, than constructs a sorted Set of them to be used as a
	 * definition for the vector space model.
	 * 
	 * @param prot
	 * @return SortedSet<String> of the respective InterproIDs in their natural
	 *         order
	 */
	public static SortedSet<String> constructVectorSpaceModel(Protein prot) {
		SortedSet<String> vectorSpaceModel = new TreeSet<String>();

		return vectorSpaceModel;
	}

	public DomainScoreCalculator(Protein protein) {
		super();
		setProtein(protein);
	}

	public Protein getProtein() {
		return protein;
	}

	public void setProtein(Protein protein) {
		this.protein = protein;
	}

}
