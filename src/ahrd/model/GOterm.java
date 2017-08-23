package ahrd.model;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

public class GOterm implements java.io.Serializable {

	// Used during deserialization to verify that the sender and receiver of a
	// serialized object have loaded classes for that object that are compatible
	// with respect to serialization.
	private static final long serialVersionUID = 1L;
	private String accession;
	private String name;
	private String ontology;
	private Boolean obsolete = false;
	private Integer annotationCount = 0;
	private Double annotationFrequency = 0.0;
	private Double informationContent = Double.POSITIVE_INFINITY; 
	/**
	 * Ancestral Gene Ontology (GO) terms of a particular GOterm include the GOterm instance itself.
	 * A GOterm is also parental to itself.
	 */
	private Set<GOterm> ancestry = new HashSet<GOterm>();

	public GOterm(String accession, String name, String ontology) {
		super();
		this.setAccession(accession);
		this.setName(name);
		this.setOntology(ontology);
	}
	
	public GOterm(String accession, String name, String ontology, Boolean obsolete) {
		this(accession, name, ontology);
		this.setObsolete(obsolete);
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

	public Boolean getObsolete() {
		return obsolete;
	}

	public void setObsolete(Boolean obsolete) {
		this.obsolete = obsolete;
	}

	public Integer getAnnotationCount() {
		return annotationCount;
	}

	public void setAnnotationCount(Integer annotationCount) {
		this.annotationCount = annotationCount;
	}

	public Double getAnnotationFrequency() {
		return annotationFrequency;
	}

	public void setAnnotationFrequency(Double annotationFrequency) {
		this.annotationFrequency = annotationFrequency;
	}

	public Double getInformationContent() {
		return informationContent;
	}

	public void setInformationContent(Double informationContent) {
		this.informationContent = informationContent;
	}

	/**
	 * Ancestral Gene Ontology (GO) terms of a particular GOterm include the GOterm instance itself.
	 * A GOterm is also parental to itself.
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
