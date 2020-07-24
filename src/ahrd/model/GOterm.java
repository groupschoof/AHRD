package ahrd.model;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * Data structure representing a more complete Gene Ontology Term. Currently not
 * in use.
 * 
 * @author Asis Hallab
 */
@Deprecated
public class GOterm {

	private String accession;
	private String name;
	private String ontology;
	/**
	 * Parental Gene Ontology (GO) term accessions include the accession of the
	 * GOterm instance. A GOterm is also parental to itself.
	 */
	private Set<String> parentAccessions;

	public GOterm(String accession, String name, String ontology) {
		super();
		this.setAccession(accession);
		this.setName(name);
		this.setOntology(ontology);
		Set<String> p = new HashSet<String>();
		p.add(accession);
		this.setParentAccessions(p);
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

	/**
	 * Parental Gene Ontology (GO) term accessions include the accession of the
	 * GOterm instance. A GOterm is also parental to itself.
	 */
	public Set<String> getParentAccessions() {
		return parentAccessions;
	}

	public void setParentAccessions(Set<String> parentAccessions) {
		this.parentAccessions = parentAccessions;
	}

}
