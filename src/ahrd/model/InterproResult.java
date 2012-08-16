package ahrd.model;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Utils.getXmlAttributeValue;
import static ahrd.controller.Utils.retrieveContentOfFirstXmlChildElement;
import static ahrd.controller.Utils.getProteinFromMemoryDatabase;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import nu.xom.Attribute;
import nu.xom.Builder;
import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Nodes;
import nu.xom.ParsingException;
import ahrd.exception.MissingInterproResultException;
import ahrd.exception.MissingProteinException;

/**
 * Representation of a Interpro-Entity of any type Domain, Family, etc.
 * 
 * @author hallab, klee
 */
public class InterproResult implements Comparable<InterproResult> {

	private String id;
	private String shortName;
	private String name;
	private String type;
	private String parentId;
	private String pfamId;
	private Set<String> contains = new HashSet<String>();
	/**
	 * The Domain-Weight is defined as the product of the two following factors.
	 * 
	 * <ul>
	 * <li>The first factor going into the similarity measurement of
	 * Domain-Architecture between two proteins. IAF is defined as :=
	 * log2(number of total proteins / number of proteins containing this
	 * domain)</li>
	 * <li>The second factor going into the similarity measurement of
	 * Domain-Architecture between two proteins. IV is defined as the number of
	 * distinct domain families adjacent to this domain.</li>
	 * </ul>
	 * 
	 * @note: see http://www.biomedcentral.com/1471-2105/10/S15/S5
	 */
	private Double domainWeight;

	private static Map<String, InterproResult> interproDb = new HashMap<String, InterproResult>();
	private static Map<String, Double> pfamDomainWeights = new HashMap<String, Double>();

	public InterproResult(String id, String shortName, String type) {
		super();
		setId(id);
		setShortName(shortName);
		setType(type);
	}

	public InterproResult(String id, String shortName, String type,
			Double domainWeight) {
		this(id, shortName, type);
		setDomainWeight(domainWeight);
	}

	public static void initialiseInterproDb() throws IOException,
			ParsingException {
		Builder parser = new Builder();
		Document doc = parser.build(new BufferedInputStream(
				new FileInputStream(new File(getSettings()
						.getPathToInterproDatabase()))));
		Nodes ipr_nodes = doc.query("//interpro");
		for (int i = 0; i < ipr_nodes.size(); i++) {
			Element ipr_el = (Element) ipr_nodes.get(i);
			InterproResult ipr = new InterproResult(getXmlAttributeValue(
					ipr_el, "id"), getXmlAttributeValue(ipr_el, "short_name"),
					getXmlAttributeValue(ipr_el, "type"));
			ipr.setName(retrieveContentOfFirstXmlChildElement(ipr_el, "name"));
			// Retrieve Parent-Id:
			Nodes resultNodes = ipr_el.query("parent_list");
			if (resultNodes.size() > 0) {
				resultNodes = resultNodes.get(0).query("rel_ref");
				if (resultNodes.size() > 0) {
					Attribute relRefAttr = ((Element) resultNodes.get(0))
							.getAttribute("ipr_ref");
					if (relRefAttr != null)
						ipr.setParentId(relRefAttr.getValue());
				}
			}
			// Initialise contained InterproResults:
			resultNodes = ipr_el.query("contains");
			if (resultNodes.size() > 0) {
				resultNodes = ((Element) resultNodes.get(0)).query("rel_ref");
				if (resultNodes.size() > 0) {
					for (int j = 0; j < resultNodes.size(); j++) {
						Element relRefElmnt = (Element) resultNodes.get(j);
						Attribute relRefAttr = relRefElmnt
								.getAttribute("ipr_ref");
						if (relRefAttr != null) {
							ipr.getContains().add(relRefAttr.getValue());
						}
					}
				}
			}
			// Add new InterproResult to the Interpro-Memory-Database:
			getInterproDb().put(ipr.getId(), ipr);
		}
	}

	/**
	 * Assigns each conserved protein domain in the memory-database its
	 * domain-weight as defined above. See mentioned article for details.
	 * Depending on AHRD's input this is done for either Pfam or InterPro
	 * domains. In the latter case this method has to be called after the memory
	 * Interpro-Database has been initialized, of course.
	 * 
	 * @throws IOException
	 * @throws NumberFormatException
	 * @Note: The input file is defined in Settings and expected to be a
	 *        tab-delimited file, in which the first column holds the
	 *        [Interpro|Pfam]-ID and the eighth column holds the domain-weight
	 *        for Eukaryotes.
	 */
	public static void parseDomainWeights() throws NumberFormatException,
			IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(
				new FileInputStream(getSettings()
						.getPathToDomainWeightsDatabase())));
		String domainId;
		String[] entry = null;
		for (String line; (line = reader.readLine()) != null;) {
			entry = line.split("\t");
			domainId = entry[0];
			Double domainWeight = Double.parseDouble(entry[7]);
			if (getSettings()
					.isDomainArchitectureSimilarityBasedOnPfamAnnotations()) {
				getPfamDomainWeights().put(domainId, domainWeight);
			} else {
				InterproResult interproEntry = getInterproDb().get(domainId);
				interproEntry.setDomainWeight(domainWeight);
			}
		}
	}

	/**
	 * Reads in a raw Interpro-Scan result file and assigns iteratively
	 * InterproResult-instances to the Proteins, specified by their
	 * Gene-Accessions.
	 * 
	 * @param proteinDb
	 * @throws IOException
	 * @throws MissingProteinException
	 */
	public static void parseInterproResult(Map<String, Protein> proteinDb)
			throws IOException, MissingProteinException {
		Set<String> missingInterproIds = new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(
				getSettings().getPathToInterproResults())));
		String iterLine = null;
		while ((iterLine = br.readLine()) != null) {
			// In case domain architecture similarities are computed based on
			// Pfam annotations parse those:
			if (getSettings().isToComputeDomainSimilarities()
					&& getSettings().getComputeDomainSimilarityOn() != null
					&& getSettings().getComputeDomainSimilarityOn().equals(
							"pfam")) {
				String pfamRegEx = "(\\S+)\\s+.*\\s(PF\\d{5})\\s.*";
				String[] annotation = parseDomainAnnotation(iterLine, pfamRegEx);
				if (annotation.length == 2) {
					getProteinFromMemoryDatabase(annotation[0], proteinDb)
							.getPfamResults().add(annotation[1]);
				}
			}
			// In order to print out conserved protein domains for the query
			// proteins, we need the InterPro annotations even in case we use
			// the Pfam annotations to compute domain architecture similarities.
			String iprRegEx = "(\\S+)\\s+.*\\s(IPR\\d{6})\\s.*";
			String[] annotation = parseDomainAnnotation(iterLine, iprRegEx);
			if (annotation.length == 2) {
				Protein prot = getProteinFromMemoryDatabase(annotation[0],
						proteinDb);
				InterproResult ipr = null;
				// WARN, if an Interpro-Result is not found in the
				// memory-database:
				if (getInterproDb().containsKey(annotation[1])) {
					ipr = getInterproDb().get(annotation[1]);
				} else {
					missingInterproIds.add(annotation[1]);
				}
				if (prot != null && ipr != null) {
					prot.getInterproResults().add(ipr);
				}
			}

		}
		if (missingInterproIds.size() > 0)
			System.err
					.println("Could not find the following Interpro-IDs in Database:\n"
							+ missingInterproIds);
	}

	/**
	 * Parses a single line from an Interpro-Scan result file extracting both
	 * the annotated protein's accession and the Interpro or Pfam identifier it
	 * has been annotated with. Which is extracted depends on the argument
	 * regex.
	 * 
	 * @Note: The argument regex should define two match groups. The first one
	 *        should match the query protein's accession and the second one the
	 *        conserved domain's identifier (Pfam or InterPro).
	 * 
	 * @param lineFromIprScanResult
	 * @return String[], in which the first element is the Protein-Accession and
	 *         the second is the Pfam-ID. Returns an <em>EMPTY</em> String
	 *         array, if line does not match expected format.
	 */
	public static String[] parseDomainAnnotation(String lineFromIprScanResult,
			String regex) {
		String[] result = new String[] {};
		Pattern pn = Pattern.compile(regex);
		Matcher ma = pn.matcher(lineFromIprScanResult);
		if (ma.matches()) {
			String protAcc = ma.group(1);
			String domainId = ma.group(2);
			if (protAcc != null && domainId != null && !protAcc.equals("")
					&& !domainId.equals("")) {
				result = new String[] { protAcc, domainId };
			}
		}
		return result;
	}

	public static Map<String, InterproResult> getInterproDb() {
		return interproDb;
	}

	public static void setInterproDb(Map<String, InterproResult> interproDb) {
		InterproResult.interproDb = interproDb;
	}

	/**
	 * Filters out all those Protein's InterproResults, who are parents of or
	 * contained by any other of the protein's InterproResults. Unfortunately
	 * this is of O(n^2)!
	 */
	public static void filterForMostInforming(Protein p)
			throws MissingInterproResultException {
		Set<InterproResult> mostInformatives = new HashSet<InterproResult>(p
				.getInterproResults());
		for (InterproResult iprToValidate : p.getInterproResults()) {
			Set<InterproResult> iprsToCompare = new HashSet<InterproResult>(p
					.getInterproResults());
			iprsToCompare.remove(iprToValidate);
			for (InterproResult iprToCompare : iprsToCompare) {
				if (iprToValidate.isParent(iprToCompare)
						|| iprToCompare.contains(iprToValidate))
					mostInformatives.remove(iprToValidate);
			}
		}
		p.setInterproResults(mostInformatives);
	}

	public int compareTo(InterproResult iprToComapre) {
		return this.getId().compareTo(iprToComapre.getId());
	}

	public boolean contains(InterproResult container)
			throws MissingInterproResultException {
		boolean isContained = false;
		if (getInterproDb().containsKey(container.getId())) {
			isContained = this.getContains().contains(container.getId());
			// Try to find container recursively
			// in this Instance's contained InterproResults:
			if (!isContained && !this.getContains().isEmpty()) {
				for (String iterContainedId : this.getContains()) {
					if (!getInterproDb().containsKey(iterContainedId))
						throw new MissingInterproResultException(
								"Could not find Interpro-Result for ID '"
										+ iterContainedId
										+ "' in Memory-Database.");
					isContained = getInterproDb().get(iterContainedId)
							.contains(container);
					if (isContained)
						return true;
				}
			}
		}
		return isContained;
	}

	/**
	 * Recursively infers if the argument InterproResult is among the ancestors
	 * of this InterproResult.
	 */
	public boolean isParent(InterproResult parent)
			throws MissingInterproResultException {
		boolean isParent = false;
		if (getInterproDb().containsKey(parent.getId())) {
			isParent = getParentId() != null
					&& getParentId().equals(parent.getId());
			// try recursive search
			if (!isParent && getParentId() != null) {
				if (!getInterproDb().containsKey(this.getParentId()))
					throw new MissingInterproResultException(
							"Could not find Interpro-Result for ID '"
									+ this.getParentId()
									+ "' in Memory-Database.");
				isParent = getInterproDb().get(this.getParentId()).isParent(
						parent);
			}
		}
		return isParent;
	}

	/**
	 * Get id.
	 * 
	 * @return id as String.
	 */
	public String getId() {
		return id;
	}

	/**
	 * Set id.
	 * 
	 * @param id
	 *            the value to set.
	 */
	public void setId(String id) {
		this.id = id;
	}

	/**
	 * Get shortName.
	 * 
	 * @return shortName as String.
	 */
	public String getShortName() {
		return shortName;
	}

	/**
	 * Set shortName.
	 * 
	 * @param shortName
	 *            the value to set.
	 */
	public void setShortName(String shortName) {
		this.shortName = shortName;
	}

	/**
	 * Get name.
	 * 
	 * @return name as String.
	 */
	public String getName() {
		return name;
	}

	/**
	 * Set name.
	 * 
	 * @param name
	 *            the value to set.
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * Get type.
	 * 
	 * @return type as String.
	 */
	public String getType() {
		return type;
	}

	/**
	 * Set type.
	 * 
	 * @param type
	 *            the value to set.
	 */
	public void setType(String type) {
		this.type = type;
	}

	/**
	 * Get parentId.
	 * 
	 * @return parentId as String.
	 */
	public String getParentId() {
		return parentId;
	}

	/**
	 * Set parentId.
	 * 
	 * @param parentId
	 *            the value to set.
	 */
	public void setParentId(String parentId) {
		this.parentId = parentId;
	}

	public String getPfamId() {
		return pfamId;
	}

	public void setPfamId(String pfamId) {
		this.pfamId = pfamId;
	}

	/**
	 * Get contains.
	 * 
	 * @return contains as Set<String>.
	 */
	public Set<String> getContains() {
		return contains;
	}

	/**
	 * Set contains.
	 * 
	 * @param contains
	 *            the value to set.
	 */
	public void setContains(Set<String> contains) {
		this.contains = contains;
	}

	public Double getDomainWeight() {
		return domainWeight;
	}

	public void setDomainWeight(Double domainWeight) {
		this.domainWeight = domainWeight;
	}

	public static Map<String, Double> getPfamDomainWeights() {
		return pfamDomainWeights;
	}

	public static void setPfamDomainWeights(
			Map<String, Double> pfamDomainWeights) {
		InterproResult.pfamDomainWeights = pfamDomainWeights;
	}

}
