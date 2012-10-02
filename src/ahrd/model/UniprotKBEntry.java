package ahrd.model;

import static ahrd.controller.Utils.retrieveAttribteValuesOfXmlChildrenElements;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.Callable;

import ahrd.exception.UniprotWebServiceAccessException;

import nu.xom.Builder;
import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Nodes;
import nu.xom.ParsingException;
import nu.xom.ValidityException;
import nu.xom.XPathContext;

public class UniprotKBEntry {

	/**
	 * In order to speed up loading of UniprotKBEntries from the RESTful service
	 * do it in parallel.
	 * 
	 * @author Asis Hallab, Mythri Bangalore
	 */
	public static class ParallelLoader implements Callable<Boolean> {

		private String accession;

		public ParallelLoader(String accession) {
			super();
			this.accession = accession;
		}

		public Boolean call() throws UniprotWebServiceAccessException {
			if (!DomainScoreCalculator.getBlastResultAccessionsToInterproIds()
					.containsKey(this.accession)
					&& !DomainScoreCalculator
							.getBlastResultAccessionsToPfamIds().containsKey(
									this.accession)) {
				String url = "NOT INITIALIZED";
				try {

					url = UniprotKBEntry.url(this.accession);
					UniprotKBEntry result = UniprotKBEntry.fromUrl(url,
							this.accession);
					DomainScoreCalculator
							.getBlastResultAccessionsToInterproIds().put(
									result.getAccession(),
									result.getIprAnnotations());
					DomainScoreCalculator.getBlastResultAccessionsToPfamIds()
							.put(result.getAccession(),
									result.getPfamAnnotations());
				} catch (Exception e) {
					System.err
							.println("Failed to access Uniprot RESTful Web Service with URL '"
									+ url
									+ "'. Probably the accession in this URL is not a Uniprot Accession?");
					throw new UniprotWebServiceAccessException(e);
				}
			}
			return true;
		}
	}

	// END of subclass!

	public static final String baseUniprotKBUrl = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/#ACCESSION#/xml";
	public static final String ACCESSION_PLACEHOLDER = "#ACCESSION#";

	public static String url(String accession)
			throws UnsupportedEncodingException {
		return baseUniprotKBUrl.replace(ACCESSION_PLACEHOLDER,
				URLEncoder.encode(accession, "UTF-8"));
	}

	/**
	 * @param url
	 * @param accession
	 *            - Provided as argument to store the accession used in the
	 *            context of this Program, which might differ from Uniprot's
	 *            latest accession.
	 * @return UniprotKBEntry as instantiated from the downloaded content.
	 * @throws IOException
	 * @throws ValidityException
	 * @throws ParsingException
	 */
	public static UniprotKBEntry fromUrl(String url, String accession)
			throws IOException, ValidityException, ParsingException {
		UniprotKBEntry result = null;
		Builder parser = new Builder();
		XPathContext c = new XPathContext("xmlns", "http://uniprot.org/uniprot");
		Document doc = parser.build(url);
		Nodes nds = doc.query("//xmlns:entry", c);
		if (nds.size() > 0) {
			Element uni = (Element) nds.get(0);
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
