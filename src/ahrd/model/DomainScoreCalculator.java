package ahrd.model;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
	 * Result from parsing SIMAP's concatenated and preprocessed feature-files.
	 * Preprocessing involves substitution of SIMAP-Hashes with the original
	 * Protein-Accessions.
	 * 
	 * @note: See awk-scripts in directory helper_scripts.
	 */
	private static Map<String, Set<String>> blastResultAccessionsToInterproIds;
	private static int dotProduct;
	private static int magnitude;

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
		String line = null;
		while ((line = reader.readLine()) != null) {
			Pattern pn = Pattern.compile("(\\S+)\\s+.*\\s(IPR\\d{6})\\s.*");
			Matcher mr = pn.matcher(line);
			if (mr.matches()) {
				accession = mr.group(1);
				interproId = mr.group(2);

				// Parse lines exactly like in method
				// InterproResult.parseInterproResult(Map<String, Protein>)
				// You should extract two values: Protein-Accession and
				// Interpro-ID
				// put those results in Map 'blastResultAccessionsToInterproIds'

				if (!getBlastResultAccessionsToInterproIds().containsKey(
						accession)) {
					getBlastResultAccessionsToInterproIds().put(accession,
							new HashSet<String>());
				}
				getBlastResultAccessionsToInterproIds().get(accession).add(
						interproId);
			}
		}
		reader.close();
	}

	private Protein protein;

	public static void constructDomainWeightVectors(Protein prot) {
		// 1.) Construct the vector space model for the Protein and its
		// BlastResults
		// 2.) Construct the domain-weights vector for the Protein itself
		// 3.) ...and all of its BlastResults
		// - Remarks: Save the above domain-weight vectors in both the Protein
		// and all BlastResults. Look at the Protein and BlastResult classes,
		// they both have a List 'domainWeights', that has to be filled with
		// instances of class Double. Do not forget to initialize those Lists
		// with new Vector() before trying to add the weights, you will generate
		// NullPointerExceptions this way. Also keep in mind, that the
		// domain-weight vectors of BlastResults, whose accession does not
		// appear in the memory database 'BlastResultAccessionsToInterproIds'
		// should be initialized to the ZERO vector, where the ZERO vector is
		// (0.0, 0.0, 0.0, ... , 0.0).

		// Map<String, Set<String>> blastToIPR =
		// getBlastResultAccessionsToInterproIds();

		List<Double> prVec = new Vector<Double>();
		List<Double> brVec = new Vector<Double>();

		SortedSet<String> vsm = constructVectorSpaceModel(prot);
		for (Iterator<String> it = vsm.iterator(); it.hasNext();) {
			String ipr = it.next();
			InterproResult interproEntry = InterproResult.getInterproDb().get(
					ipr);
			double weight = interproEntry.getDomainWeight();

			if (prot.getInterproResults().contains(interproEntry)) {
				prVec.add(weight);
			} else
				prVec.add(0.0);
		}
		prot.setDomainWeights(prVec);

		for (String blastDb : prot.getBlastResults().keySet()) {
			for (BlastResult br : prot.getBlastResults().get(blastDb)) {
				Set<String> iprSet = getBlastResultAccessionsToInterproIds()
						.get(br.getAccession());

				if (iprSet != null && iprSet.size() > 0) {
					for (Iterator<String> it2 = iprSet.iterator(); it2
							.hasNext();) {
						String ipr = it2.next();
						InterproResult interproEntry = InterproResult
								.getInterproDb().get(ipr);
						double weight = interproEntry.getDomainWeight();

						if (vsm.contains(interproEntry.getId())) {
							brVec.add(weight);
						} else
							brVec.add(0.0);
					}
				}
				if (!getBlastResultAccessionsToInterproIds().containsKey(
						br.getAccession())) {
					brVec = new Vector<Double>(Arrays.asList(0.0, 0.0, 0.0,
							0.0, 0.0, 0.0));
				}
				br.setDomainWeights(brVec);
			}
		}
	}

	/**
	 * Calculates the cosine of angle between the two argument vectors as a
	 * measure of their similarity: sim(x,y) = dot-product(x,y) / (||x||*||y||)
	 * 
	 * @param x
	 * @param y
	 * @return sim(x,y)
	 */

	public static Double domainWeightSimilarity(List<Double> prVec,
			List<Double> brVec) {

		// According to the above mentioned article, calculate the cosine of the
		// angle between between the two argument vectors, using the dot-product
		// in the numerator and the product of euclidean lengths as the
		// denominator.
		Double dotProduct = 0.0;
		for (int i = 0; i < prVec.size(); i++) {
			dotProduct += (prVec.get(i) * brVec.get(i));
		}

		Double magPr = 0.0;
		Double magBr = 0.0;
		for (int i = 0; i < prVec.size(); i++) {
			magPr += Math.pow(prVec.get(i), 2);
			magBr += Math.pow(brVec.get(i), 2);
		}
		Double magnitude = Math.sqrt(magPr) * Math.sqrt(magBr);

		Double dws = dotProduct / magnitude;
		return dws;
	}

	/**
	 * 
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
		// 1.) Construct the template for the Vector of domain-weights!
		for (InterproResult ir : prot.getInterproResults()) {
			vectorSpaceModel.add(ir.getId());
		}
		// For each BlastResult, check, if there's an entry in the memory
		// database 'BlastResultAccessionsToInterproIds' using the BlastResult's
		// Accession as lookup key. Add all found assigned Interpro-Ids to the
		// vector space model.
		for (String blastDb : prot.getBlastResults().keySet()) {
			for (BlastResult br : prot.getBlastResults().get(blastDb)) {
				br.getAccession();
				if (getBlastResultAccessionsToInterproIds().containsKey(
						br.getAccession())) {
					vectorSpaceModel.addAll(blastResultAccessionsToInterproIds
							.get(br.getAccession()));
				}
			}

		}
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
