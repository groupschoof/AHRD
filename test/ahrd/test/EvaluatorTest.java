package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.sql.SQLException;

import org.junit.Before;
import org.junit.Test;
import org.xml.sax.SAXException;

import ahrd.controller.Evaluator;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingInterproResultException;
import ahrd.exception.MissingProteinException;
import ahrd.model.GOterm;

public class EvaluatorTest {

	private Evaluator evaluator;

	@Before
	public void setup() throws IOException {
		this.evaluator = new Evaluator("./test/resources/evaluator_example_hrd_and_go.yml");
	}

	@Test
	public void testSetupGroundTruthDescriptions() throws IOException,
			MissingAccessionException {
		evaluator.initializeProteins();
		evaluator.setupGroundTruthDescriptions();
		assertEquals(3, evaluator.getProteins().size());
		assertNotNull(evaluator.getProteins()
				.get("gene:chr01.1056:mRNA:chr01.1056")
				.getEvaluationScoreCalculator().getGroundTruthDescription());
		assertNotNull(evaluator.getProteins()
				.get("gene:chr01.502:mRNA:chr01.502")
				.getEvaluationScoreCalculator().getGroundTruthDescription());
		assertEquals("Receptor-like protein kinase", evaluator.getProteins()
				.get("gene:chr01.1056:mRNA:chr01.1056")
				.getEvaluationScoreCalculator().getGroundTruthDescription()
				.getDescription());
		assertEquals(4,
				evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056")
						.getEvaluationScoreCalculator()
						.getGroundTruthDescription().getTokens().size());
	}

	@Test
	public void testTokenizesUnchangedBlastResults() throws IOException,
			MissingAccessionException, MissingProteinException, SAXException {
		evaluator.initializeProteins();
		evaluator.parseBlastResults();
		// TEST:
		assertEquals(
				"Parsing the BlastResults should also trigger tokenization of the Description-Lines of the unchanged BlastHits!",
				4, evaluator.getProteins().get("gene:chr01.502:mRNA:chr01.502")
						.getEvaluationScoreCalculator()
						.getBestUnchangedBlastResults().get("swissprot")
						.getTokens().size());
	}
	
	@Test
	public void testSetupCompetitors() throws IOException, MissingAccessionException, MissingInterproResultException, SQLException {
		evaluator.initializeProteins();
		evaluator.setupGroundTruthDescriptions();
		evaluator.assignHumanReadableDescriptions();
		evaluator.setupGoAnnotationEvaluation();
		evaluator.setupCompetitors();
		
		assertEquals(3, evaluator.getProteins().size());
		
		assertNotNull(
				"After setting up Blast2GoAnnots the Evaluator should have assigned a Blast2GoAnnot to the protein with accession 'gene:chr01.502:mRNA:chr01.502'.",
				evaluator.getProteins().get("gene:chr01.502:mRNA:chr01.502")
						.getEvaluationScoreCalculator().getCompetitorAnnotations());
		
		assertNotNull(
				"After setting up Competitors the Evaluator should have assigned a CompetitorAnnotation to the protein with accession 'gene:chr01.1056:mRNA:chr01.1056'.",
				evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056")
						.getEvaluationScoreCalculator().getCompetitorAnnotations());
		assertNotNull(evaluator.getProteins()
				.get("gene:chr01.1056:mRNA:chr01.1056")
				.getEvaluationScoreCalculator().getCompetitorAnnotations());
		assertEquals(1, evaluator.getProteins()
				.get("gene:chr01.1056:mRNA:chr01.1056")
				.getEvaluationScoreCalculator().getCompetitorAnnotations().size());
		assertEquals("nrpb6a dna binding dna-directed rna polymerase",
				evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056")
						.getEvaluationScoreCalculator().getCompetitorAnnotations().get("blast2go").getDescription());
		assertEquals(6,
				evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056")
						.getEvaluationScoreCalculator().getCompetitorAnnotations().get("blast2go").getEvaluationTokens().size());
		assertEquals(2, evaluator.getProteins()
				.get("gene:chr01.1056:mRNA:chr01.1056")
				.getEvaluationScoreCalculator().getCompetitorAnnotations().get("blast2go").getGoAnnotations().size());
		for (GOterm term : evaluator.getProteins()
				.get("gene:chr01.1056:mRNA:chr01.1056")
				.getEvaluationScoreCalculator().getCompetitorAnnotations().get("blast2go").getGoAnnotations()) {
			assertTrue(term.getAccession().equals("GO:0005525") || term.getAccession().equals("GO:0007264"));
		}
		
		assertNotNull(evaluator.getProteins()
				.get("P04637")
				.getEvaluationScoreCalculator().getCompetitorAnnotations());
		assertEquals(1, evaluator.getProteins()
				.get("P04637")
				.getEvaluationScoreCalculator().getCompetitorAnnotations().size());
		assertEquals("", evaluator.getProteins()
				.get("P04637")
				.getEvaluationScoreCalculator().getCompetitorAnnotations().get("blast2go").getDescription());
		assertEquals(3, evaluator.getProteins()
				.get("P04637")
				.getEvaluationScoreCalculator().getCompetitorAnnotations().get("blast2go").getGoAnnotations().size());
		for (GOterm term : evaluator.getProteins()
				.get("P04637")
				.getEvaluationScoreCalculator().getCompetitorAnnotations().get("blast2go").getGoAnnotations()) {
			assertTrue(term.getAccession().equals("GO:0005524") || term.getAccession().equals("GO:0006284") || term.getAccession().equals("GO:0000733"));
		}
	}
}
