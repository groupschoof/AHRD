package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.sql.SQLException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import nu.xom.ParsingException;

import org.junit.Before;
import org.junit.Test;
import org.semanticweb.owlapi.model.OWLOntologyCreationException;
import org.xml.sax.SAXException;

import ahrd.controller.AHRD;
import ahrd.controller.Utils;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;
import ahrd.model.Protein;
import ahrd.model.ReferenceGoAnnotation;

public class GoAnnotationReferenceTest {

	private AHRD ahrd;

	@Before
	public void setUp() throws IOException {
		ahrd = new AHRD(
				"./test/resources/ahrd_input_seq_sim_table_go_prediction.yml");
	}

	@Test
	public void testHasGeneOntologyAnnotations() {
		// Should have GO-Annotations with default test-Settings:
		assertTrue(getSettings().hasGeneOntologyAnnotations());
		getSettings().removeAllPathToGeneOntologyReferemces();
		assertTrue(!getSettings().hasGeneOntologyAnnotations());
		getSettings().setPathToGeneOntologyReference("swissprot","/not/existing/path.raw");
		assertTrue(!getSettings().hasGeneOntologyAnnotations());
		getSettings().setPathToGeneOntologyReference("swissprot","./test/resources/reference_gene_ontology_annotations_uniprotKB_GOA.txt");
		assertTrue(getSettings().hasGeneOntologyAnnotations());
	}

	@Test
	public void testUniqueShortAccessions() throws IOException,
			MissingAccessionException, MissingProteinException, SAXException,
			ParsingException {
		ahrd.setup(false);
		assertNotNull(ahrd.getUniqueBlastResultShortAccessions());
		// Somehow assertEquals does not work on Collections as expected, hence
		// the following work-around:
		List<String> expectedUniqueShortAccessions = Utils
				.fromFile("./test/resources/all_blast_hits_blast8_tabular_searches.txt");
		for (String sa : expectedUniqueShortAccessions) {
			assertTrue(ahrd.getUniqueBlastResultShortAccessions().contains(sa));
		}
		assertEquals(expectedUniqueShortAccessions.size(), ahrd
				.getUniqueBlastResultShortAccessions().size());
		ahrd.getUniqueBlastResultShortAccessions().removeAll(
				expectedUniqueShortAccessions);
		assertEquals(0, ahrd.getUniqueBlastResultShortAccessions().size());
	}

	@Test
	public void testParseGoAnnotationReference() throws IOException,
			MissingAccessionException, MissingProteinException, SAXException,
			ParsingException {
		ahrd.setup(false);
		assertNotNull(ahrd.getGoAnnotationReference());
		assertTrue(!ahrd.getGoAnnotationReference().isEmpty());
		assertEquals(4, ahrd.getGoAnnotationReference().size());
		Set<ReferenceGoAnnotation> reference = ahrd.getGoAnnotationReference().get("AT1G01040");
		Set<String> dbGos = new HashSet<String>();
		for (ReferenceGoAnnotation annotation : reference) {
			dbGos.add(annotation.getGoTerm());
		}
		assertTrue(dbGos.contains("GO:0005634"));
		assertTrue(dbGos.contains("GO:0008026"));
	}

	@Test
	public void testAnnotatesGoTerms() throws IOException,
			MissingAccessionException, MissingProteinException, SAXException,
			ParsingException, SQLException, OWLOntologyCreationException {
		ahrd.setup(false);
		ahrd.assignGeneOntologyTerms();
		Protein p = ahrd.getProteins().get("gene:chr01.1056:mRNA:chr01.1056");
		assertNotNull(p.getGoResults());
		assertEquals(2, p.getGoResults().size());
		assertTrue(p.getGoResults().contains("GO:0006355"));
		assertTrue(p.getGoResults().contains("GO:0043401"));
	}
}
