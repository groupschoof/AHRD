package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.junit.Before;
import org.junit.Test;

import ahrd.controller.AHRD;
import ahrd.model.InterproResult;

public class RunAhrdWithDomainArchitectureTest {

	private AHRD ahrd;

	@Before
	public void setUp() throws IOException {
		ahrd = new AHRD("./test/resources/ahrd_dom_arch_sim_input.yml");
	}

	/**
	 * Check each step of AHRD executed with options set to take domain
	 * architecture similarity (DAS) into account. DAS is computed on InterPro
	 * domains.
	 * 
	 * @throws Exception
	 */
	@Test
	public void testRunInterpro() throws Exception {
		// Check correct initialization:
		ahrd.setup(false);

		// Verify the expected InterPro domains have been loaded:
		assertNotNull(
				"InterPro Entry 'IPR000001' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR000001"));
		assertNotNull(
				"InterPro Entry 'IPR000003' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR000003"));
		assertNotNull(
				"InterPro Entry 'IPR000006' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR000006"));
		assertNotNull(
				"InterPro Entry 'IPR000535' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR000535"));
		assertNotNull(
				"InterPro Entry 'IPR000536' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR000536"));
		assertNotNull(
				"InterPro Entry 'IPR001723' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR001723"));
		assertNotNull(
				"InterPro Entry 'IPR001732' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR001732"));
		assertNotNull(
				"InterPro Entry 'IPR001733' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR001733"));
		assertNotNull(
				"InterPro Entry 'IPR003019' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR003019"));
		assertNotNull(
				"InterPro Entry 'IPR008946' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR008946"));
		assertNotNull(
				"InterPro Entry 'IPR013806' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR013806"));
		assertNotNull(
				"InterPro Entry 'IPR016040' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR016040"));
		assertNotNull(
				"InterPro Entry 'IPR019756' was not correctly loaded from the interpro database!",
				InterproResult.getInterproDb().get("IPR019756"));

		// Domain Weights Database has to have been parsed:
		InterproResult ipr = InterproResult.getInterproDb().get("IPR000001");
		Double dw = ipr.getDomainWeight();
		assertNotNull(
				"InterPro Entry 'IPR000001' has no Domain Weight assigned!", dw);
		assertEquals(0.189433136086726, dw, 0.0);
	}

}
