package ahrd.model;

import static ahrd.controller.Utils.retrieveAttribteValuesOfXmlChildrenElements;
import static ahrd.controller.Utils.retrieveContentOfFirstXmlChildElement;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import nu.xom.Builder;
import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Nodes;
import nu.xom.ParsingException;
import nu.xom.ValidityException;
import nu.xom.XPathContext;

public class UniprotKBEntry {

	/**
	 * In order to speed up loading of UniprotKBEntries from the RESTful
	 * service, do it in parallel.
	 * 
	 * @author Asis Hallab
	 */
	public class ParallelLoader implements Runnable {

		private String accession;

		public ParallelLoader(String accession) {
			super();
			this.accession = accession;
		}

		public void run() {
			if (!DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
					.containsKey(this.accession)) {
				try {
					UniprotKBEntry result = UniprotKBEntry
							.fromUrl(this.accession);
					// todo

				} catch (Exception e) {
					e.printStackTrace(System.err);
				}
			}
		}
	}

	// END of subclass!

	public static final String baseUniprotKBUrl = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/#ACCESSION#/xml";
	public static final String ACCESSION_PLACEHOLDER = "#ACCESSION#";

	public static String url(String accession) {
		return baseUniprotKBUrl.replace(ACCESSION_PLACEHOLDER, accession);
	}

	public static UniprotKBEntry fromUrl(String url) throws IOException,
			ValidityException, ParsingException {
		UniprotKBEntry result = null;
		Builder parser = new Builder();
		XPathContext c = new XPathContext("xmlns", "http://uniprot.org/uniprot");
		Document doc = parser.build(url);
		Nodes nds = doc.query("//xmlns:entry", c);
		if (nds.size() > 0) {
			Element uni = (Element) nds.get(0);
			String accession = retrieveContentOfFirstXmlChildElement(uni,
					"xmlns:accession", c);
			// Instantiate new UniprotKBEntry, if and only if we find a valid
			// accession:
			if (accession != null && !accession.equals("")) {
				result = new UniprotKBEntry(accession);
				// Add Interpro-Annotations
				result.setIprAnnotations(retrieveAttribteValuesOfXmlChildrenElements(
						uni, "xmlns:dbReference[@type='InterPro']", "id", c));
				// Add PFAM-Annotations
				result.setPfamAnnotations(retrieveAttribteValuesOfXmlChildrenElements(
						uni, "xmlns:dbReference[@type='Pfam']", "id", c));
			}
		}
		return result;
	}

	private String accession;
	private Set<String> iprAnnotations = new HashSet<String>();
	private Set<String> pfamAnnotations = new HashSet<String>();

	public UniprotKBEntry() {
		super();
	}

	public UniprotKBEntry(String accession) {
		super();
		setAccession(accession);
	}

	public String getAccession() {
		return accession;
	}

	public void setAccession(String accession) {
		this.accession = accession;
	}

	public Set<String> getIprAnnotations() {
		return iprAnnotations;
	}

	public void setIprAnnotations(Set<String> iprAnnotations) {
		this.iprAnnotations = iprAnnotations;
	}

	public Set<String> getPfamAnnotations() {
		return pfamAnnotations;
	}

	public void setPfamAnnotations(Set<String> pfamAnnotations) {
		this.pfamAnnotations = pfamAnnotations;
	}

}
