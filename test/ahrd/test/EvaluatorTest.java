package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static ahrd.controller.DatabaseSetup.createOrUpdateAhrdDatabase;

import java.io.IOException;

import org.junit.Before;
import org.junit.Test;
import org.xml.sax.SAXException;

import ahrd.controller.Evaluator;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;
import ahrd.model.Blast2GoAnnot;

public class EvaluatorTest {

	private Evaluator evaluator;

	@Before
	public void setup() throws IOException {
		this.evaluator = new Evaluator("./test/resources/ahrd_input.yml");
	}

	@Test
	public void testSetupBlast2GoAnnots() throws IOException, MissingAccessionException {
		evaluator.initializeProteins();
		evaluator.setupBlast2GoAnnots();
		assertEquals(2, evaluator.getProteins().size());
		assertNotNull(
				"After setting up Blast2GoAnnots the Evaluator should have assigned a Blast2GoAnnot to the protein with accession 'gene:chr01.1056:mRNA:chr01.1056'.",
				evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056").getEvaluationScoreCalculator()
						.getBlast2GoAnnots());
		assertNotNull(
				"After setting up Blast2GoAnnots the Evaluator should have assigned a Blast2GoAnnot to the protein with accession 'gene:chr01.1056:mRNA:chr01.1056'.",
				evaluator.getProteins().get("gene:chr01.502:mRNA:chr01.502").getEvaluationScoreCalculator()
						.getBlast2GoAnnots());
		assertEquals("nrpb6a dna binding dna-directed rna polymerase",
				evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056").getEvaluationScoreCalculator()
						.getBlast2GoAnnots().toArray(new Blast2GoAnnot[] {})[0].getDescription());
		assertEquals(6, evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056").getEvaluationScoreCalculator()
				.getBlast2GoAnnots().toArray(new Blast2GoAnnot[] {})[0].getEvaluationTokens().size());
	}

	@Test
	public void testSetupReferences() throws IOException, MissingAccessionException {
		evaluator.initializeProteins();
		evaluator.setupReferences();
		assertEquals(2, evaluator.getProteins().size());
		assertNotNull(evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056").getEvaluationScoreCalculator()
				.getReferenceDescription());
		assertNotNull(evaluator.getProteins().get("gene:chr01.502:mRNA:chr01.502").getEvaluationScoreCalculator()
				.getReferenceDescription());
		assertEquals("Receptor-like protein kinase", evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056")
				.getEvaluationScoreCalculator().getReferenceDescription().getDescription());
		assertEquals(4, evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056").getEvaluationScoreCalculator()
				.getReferenceDescription().getTokens().size());
	}

	@Test
	public void testTokenizesUnchangedBlastResults()
			throws IOException, MissingAccessionException, MissingProteinException, SAXException {
		createOrUpdateAhrdDatabase(false);
		evaluator.initializeProteins();
		evaluator.parseBlastResults();
		// TEST:
		assertEquals(
				"Parsing the BlastResults should also trigger tokenization of the Description-Lines of the unchanged BlastHits!",
				4, evaluator.getProteins().get("gene:chr01.502:mRNA:chr01.502").getEvaluationScoreCalculator()
						.getUnchangedBlastResults().get("swissprot").getTokens().size());
	}
}
