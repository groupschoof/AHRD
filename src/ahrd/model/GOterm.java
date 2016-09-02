package ahrd.model;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

public class GOterm {

	private String accession;
	private String name;
	private String ontology;
	private Integer frequency;
	private Double probability;
	private Double informationContent; 
	/**
	 * Parental Gene Ontology (GO) term accessions include the accession of the
	 * GOterm instance. A GOterm is also parental to itself.
	 */
	private Set<GOterm> ancestry;

	public GOterm(String accession, String name, String ontology) {
		super();
		this.setAccession(accession);
		this.setName(name);
		this.setOntology(ontology);
		this.setFrequency(0);
		this.setProbability(0.0);
		this.setInformationContent(Double.POSITIVE_INFINITY);
		this.ancestry = new HashSet<GOterm>();
	}

	/**
	 * Extracts all pairwise distinct Gene Ontology (GO) term accessions present
	 * in the Collection of GOterm instances in argument gts.
	 * 
	 * @param gts
	 * @return Set<String>
	 */
	public static Set<String> uniqueAccessions(Collection<GOterm> gts) {
		Set<String> res = new HashSet<String>();
		for (GOterm g : gts)
			res.add(g.getAccession());
		return res;
	}

	public String getAccession() {
		return accession;
	}

	public void setAccession(String accession) {
		this.accession = accession;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getOntology() {
		return ontology;
	}

	public void setOntology(String ontology) {
		this.ontology = ontology;
	}

	public Integer getFrequency() {
		return frequency;
	}

	public void setFrequency(Integer frequency) {
		this.frequency = frequency;
	}

	public Double getProbability() {
		return probability;
	}

	public void setProbability(Double probability) {
		this.probability = probability;
	}

	public Double getInformationContent() {
		return informationContent;
	}

	public void setInformationContent(Double informationContent) {
		this.informationContent = informationContent;
	}

	/**
	 * Parental Gene Ontology (GO) term accessions include the accession of the
	 * GOterm instance. A GOterm is also parental to itself.
	 */
	public Set<GOterm> getAncestry() {
		return ancestry;
	}

	public void setAncestry(Set<GOterm> ancestry) {
		this.ancestry = ancestry;
	}
	
	public void addTermToAncestry(GOterm term) {
		this.ancestry.add(term);
	}

}
