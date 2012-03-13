package ahrd.test;

import java.io.IOException;
import java.util.Map;

import static junit.framework.Assert.*;

import nu.xom.ParsingException;

import org.junit.Before;
import org.junit.Test;

import ahrd.model.InterproResult;
import ahrd.model.Protein;

import ahrd.exception.MissingProteinException;
import ahrd.exception.MissingInterproResultException;

public class InterproResultTest {

	public InterproResultTest() {
		super();
	}

	@Before
	public void initialiseInterproDb() throws IOException, ParsingException {
		TestUtils.initTestSettings();
		InterproResult.initialiseInterproDb();
	}

	@Test
	public void testInitialisationOfInterproDb() {
		assertNotNull(InterproResult.getInterproDb());
		assertTrue(InterproResult.getInterproDb().size() > 0);
		InterproResult ipr = InterproResult.getInterproDb().get("IPR000003");
		assertNotNull(ipr);
		assertEquals("IPR000003", ipr.getId());
		assertEquals("IPR001723", ipr.getParentId());
		assertEquals("Retinoid-X_rcpt", ipr.getShortName());
		assertEquals("Kringle", InterproResult.getInterproDb().get("IPR000001")
				.getShortName());
		assertEquals("Family", ipr.getType());
		assertEquals("Retinoid X receptor", ipr.getName());
		assertTrue(InterproResult.getInterproDb().containsKey("IPR000535"));
		assertTrue(InterproResult.getInterproDb().containsKey("IPR000536"));
	}

	@Test
	public void testRecursiveParentSearch()
			throws MissingInterproResultException {
		InterproResult irpChild = InterproResult.getInterproDb().get(
				"IPR000003");
		InterproResult irpGrandParent = InterproResult.getInterproDb().get(
				"IPR013806");
		assertTrue(irpChild.isParent(irpGrandParent));
	}

	@Test
	public void testRecursiveContainsSearch()
			throws MissingInterproResultException {
		InterproResult iprContainer = InterproResult.getInterproDb().get(
				"IPR000003");
		InterproResult iprContainee = InterproResult.getInterproDb().get(
				"IPR000535");
		InterproResult ipr1stLvlConatinee = InterproResult.getInterproDb().get(
				"IPR000536");
		assertTrue(iprContainer.contains(ipr1stLvlConatinee));
		assertTrue(iprContainer.contains(iprContainee));
	}

	@Test
	public void testParseInterproResults() throws MissingProteinException,
			MissingInterproResultException {
		Map<String, Protein> proteinDb = TestUtils.mockProteinDb();
		try {
			InterproResult.parseInterproResult(proteinDb);
		} catch (IOException e) {
			e.printStackTrace(System.out);
		}
		assertEquals(1, proteinDb.get("gene:chr01.502:mRNA:chr01.502")
				.getInterproResults().size());
		assertEquals(2, proteinDb.get("gene:chr01.1056:mRNA:chr01.1056")
				.getInterproResults().size());
		assertTrue(proteinDb.get("gene:chr01.502:mRNA:chr01.502")
				.getInterproResults().contains(
						InterproResult.getInterproDb().get("IPR000535")));
		assertTrue(proteinDb.get("gene:chr01.1056:mRNA:chr01.1056")
				.getInterproResults().contains(
						InterproResult.getInterproDb().get("IPR000006")));
		assertTrue(proteinDb.get("gene:chr01.1056:mRNA:chr01.1056")
				.getInterproResults().contains(
						InterproResult.getInterproDb().get("IPR000536")));
	}

	@Test
	public void testFilterInterproResults() throws Exception {
		Protein p = TestUtils.mockProtein();
		InterproResult ipr1 = new InterproResult("IPR:000001", "short name 1",
				"domain");
		InterproResult ipr2 = new InterproResult("IPR:000002", "short name 2",
				"domain");
		InterproResult ipr3 = new InterproResult("IPR:000003", "short name 3",
				"domain");
		InterproResult ipr4 = new InterproResult("IPR:000004", "short name 4",
				"domain");
		InterproResult ipr5 = new InterproResult("IPR:000005", "short name 5",
				"domain");
		InterproResult ipr6 = new InterproResult("IPR:000006", "short name 6",
				"domain");
		InterproResult ipr7 = new InterproResult("IPR:000007", "short name 7",
				"domain");
		// Parent -> Child:
		// ipr1 -> ipr3 -> ipr4
		ipr3.setParentId(ipr1.getId());
		ipr4.setParentId(ipr3.getId());
		// Container -> Contained:
		// ipr2 -> ipr5 + ipr6 -> ipr7
		ipr2.getContains().add(ipr5.getId());
		ipr2.getContains().add(ipr6.getId());
		ipr6.getContains().add(ipr7.getId());
		// Add InterproResults to Protein and InterproResult-Database:
		p.getInterproResults().add(ipr1);
		InterproResult.getInterproDb().put(ipr1.getId(), ipr1);
		p.getInterproResults().add(ipr2);
		InterproResult.getInterproDb().put(ipr2.getId(), ipr2);
		p.getInterproResults().add(ipr3);
		InterproResult.getInterproDb().put(ipr3.getId(), ipr3);
		p.getInterproResults().add(ipr4);
		InterproResult.getInterproDb().put(ipr4.getId(), ipr4);
		p.getInterproResults().add(ipr5);
		InterproResult.getInterproDb().put(ipr5.getId(), ipr5);
		p.getInterproResults().add(ipr6);
		InterproResult.getInterproDb().put(ipr6.getId(), ipr6);
		p.getInterproResults().add(ipr7);
		InterproResult.getInterproDb().put(ipr7.getId(), ipr7);

		// After filtering Protein p should only have the
		// InterproResults ipr1 and ipr2:
		InterproResult.filterForMostInforming(p);
		assertEquals(2, p.getInterproResults().size());
		assertTrue(p.getInterproResults().contains(ipr1));
		assertTrue(p.getInterproResults().contains(ipr2));
	}
}
