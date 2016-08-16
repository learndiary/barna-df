package barna.geneid;

import barna.commons.Execute;
import barna.geneid.GeneID;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/8/12
 * Time: 6:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProfileTest {

    @BeforeClass
    public static void initExecuter() {
        Execute.initialize(2);
    }

    @AfterClass
    public static void shutdownExecuter() {
        Execute.shutdown();
    }

    @Test
    public void testReadProfile() throws Exception {
        GeneIDsettings settings= new GeneIDsettings();
        GParam[] isochores= Profile.readParam(null, settings);
    }

    @Test
    public void testScoreSpliceSites() throws Exception {
        GeneIDsettings settings= new GeneIDsettings();
        GParam[] isochores= Profile.readParam(null, settings);
        System.currentTimeMillis();

        // sequence length: (dimesion+ order)

        // prefix_donor = (offset+ order)
        // suffix_donor = (dimension- offset- order- 2)
        // DonorProfile: order= 1, offset= 1, dimension= 9
        // prefix 2, suffix
        float donScore= GeneID.scoreDonor("GCGTACCCC", isochores[0].DonorProfile);
        System.err.println("score "+ donScore);

        // prefix_acceptor = (dimension- offset- 2)+ order
        // suffix_acceptor = offset
        // AcceptorProfile: order= 1, offset= 24, dimension= 27
        // "CTCTCTCTCTCTCTCTCTCTCTAGCGC"
        // "CTCTCTCTCTCTCTCTCTCTCTCAGCGG" => this sequence had 28 instead of 27 chars, removed 5' most C
        float accScore= GeneID.scoreAcceptor("TCTCTCTCTCTCTCTCTCTCTCAGCGG", isochores[0].AcceptorProfile, null, null);
        System.err.println("score "+ accScore);

//        myGID.buildDonors(
//                "",
//                (short) 0,
//                "", // type
//                "", // sub-type
//                isochores[0].DonorProfile, // Profile
//                new Site[1], //Site[] st
//                0,  // l1
//                0,  // l2
//                0,  // ns
//                1,  // nsites
//                GeneIDconstants.FORWARD,    // Strand
//                null    // Pack external
//        );

    }

}
