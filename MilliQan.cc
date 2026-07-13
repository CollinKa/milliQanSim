#include "mqDetectorConstruction.hh"
#include "mqPrimaryGeneratorAction.hh"
#include "MilliQPrimaryGeneratorAction.hh"
#include "mqRunAction.hh"
#include "mqEventAction.hh"
#include "mqSteppingAction.hh"
#include "mqSteppingVerbose.hh"
#include "mqTrackingAction.hh"
#include "mqStackingAction.hh"
#include "mqSession.hh"
#include "mqHistoManager.hh"

#include "mqShieldingPhysicsList.hh"
//#include "MilliQPhysicsList.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>



//#include "QGSP_BERT_HP.hh"

/*Grafikausgabe----------------------------------------------------------------*/
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

int main(int argc, char** argv) {
    boost::property_tree::ptree pt;
    boost::property_tree::ptree particlePTree;
    G4UIExecutive *uiExec;
    if (argc == 1 || argc <= 3) {
    uiExec = new G4UIExecutive(argc, argv);
    }

    // Default configuration file
    std::string configFile = "config/onepc.ini";

    // Handle arguments
    if (argc > 1) {
        configFile = argv[1];
        G4cout << "MilliQan> Using configuration file: " << configFile << G4endl;
    } else {
        G4cout << "MilliQan> Using default configuration file: " << configFile << G4endl;
    }

    // Load the main configuration file
    try {
        boost::property_tree::ini_parser::read_ini(configFile, pt);
    } catch (boost::property_tree::ptree_error& e) {
        G4ExceptionDescription msg;
        msg << G4endl << "Configuration file error: " << e.what() << G4endl;
        G4Exception("MilliQan::main()", "MilliQan::ConfigFileReadError", FatalException, msg);
    }

    // Load the particle configuration file
    try {
        boost::property_tree::ini_parser::read_ini(pt.get<std::string>("Configuration.ParticleConfigFile"), particlePTree);
    } catch (boost::property_tree::ptree_error& e) {
        G4ExceptionDescription msg;
        msg << G4endl << "Particle configuration file error: " << e.what() << G4endl;
        G4Exception("MilliQan::main()", "MilliQan::ConfigFileReadError", FatalException, msg);
    }

    // Retrieve particle properties
    G4int jobNumber = particlePTree.get<G4int>("ParticleProperties.JobNumber");
    G4int eventsPerJob = particlePTree.get<G4int>("ParticleProperties.EventsPerJob");
    G4double eventWeight = particlePTree.get<G4double>("ParticleProperties.EventWeight");
    G4int processID = particlePTree.get<G4int>("ParticleProperties.ProcessID");

    G4int eventOffset = jobNumber * eventsPerJob;

    // Initialize the randomizer for primary particles
    CLHEP::RanluxEngine defaultEngine(1234567, 4);
    G4Random::setTheEngine(&defaultEngine);

    // Set verbose output class
    G4VSteppingVerbose* verbosity = new mqSteppingVerbose;
    G4VSteppingVerbose::SetInstance(verbosity);

    // Initialize the Run Manager
    G4RunManager* runManager = new G4RunManager;

    // Set User Initialization classes
    mqDetectorConstruction* detector = new mqDetectorConstruction;
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new mqShieldingList(0, "HP", pt));

    // Visualization
#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif

    // Set Analysis Manager
    mqHistoManager* histo = new mqHistoManager();

    // Set User Action classes
    runManager->SetUserAction(new mqRunAction(histo));
    runManager->SetUserAction(new mqEventAction(histo, eventOffset, eventWeight, processID));
    runManager->SetUserAction(new mqSteppingAction(histo));
    runManager->SetUserAction(new MilliQPrimaryGeneratorAction(pt, eventOffset));
    runManager->SetUserAction(new mqTrackingAction(histo));
    runManager->SetUserAction(new mqStackingAction);

    // Get the UI manager
    G4UImanager* UI = G4UImanager::GetUIpointer();
    mqSession* LoggedSession = new mqSession();
    UI->SetCoutDestination(LoggedSession);

    // Handle different execution modes based on arguments
    /*
if (argc == 1) {
    // Error: No arguments provided
    G4ExceptionDescription msg;
    msg << "No configuration or macro file specified. Usage:" << G4endl;
    msg << "./MilliQan [configFile] [macroFile]" << G4endl;
    G4Exception("MilliQan::main()", "MilliQan::InvalidArgumentCount", FatalException, msg);
} else if (argc == 2) {
    // Configuration file provided, but no macro file
    G4cout << "MilliQan> Configuration file specified: " << argv[1] << G4endl;
    G4cout << "No macro file specified. Exiting." << G4endl;
    return 0;
    */
if (argc == 1) {
		G4cout << "MilliQan> Enter interactive mode." << G4endl;
		// Define (G)UI terminal for interactive mode
		// G4UIterminal is a (dumb) terminal
		G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
		session = new G4UIterminal(new G4UItcsh);
#else
		session = new G4UIterminal();
#endif
		// 		UI->ApplyCommand("/control/execute run25.mac");
		session->SessionStart();
		delete session;
} else if (argc <= 3) {
    // Configuration file and macro file provided
    G4cout << "MilliQan> Using configuration file: " << argv[1] << G4endl;
    G4cout << "MilliQan> Running macro file: " << argv[2] << G4endl;

    // Execute the macro file
    G4String command = "/control/execute ";
    G4String macroFile = argv[2];
    UI->ApplyCommand(command + macroFile);
    G4UIExecutive* uiExec = new G4UIExecutive(argc, argv);
    //LoggedSession->SessionStart();
    uiExec->SessionStart();

    //G4UIExecutive* uiExec = new G4UIExecutive(argc, argv);
    uiExec->SessionStart();
    delete uiExec;

    G4cout << "MilliQan> Simulation complete. Exiting." << G4endl;
    //return 0;
} else {
    // Invalid number of arguments
    G4ExceptionDescription msg;
    msg << "Invalid number of arguments. Usage:" << G4endl;
    msg << "./MilliQan [configFile] [macroFile]" << G4endl;
    G4Exception("MilliQan::main()", "MilliQan::InvalidArgumentCount", FatalException, msg);
}



    // Cleanup
#ifdef G4VIS_USE
    delete visManager;
#endif
    delete runManager;
    delete verbosity;
    delete LoggedSession;
    delete histo;

    return 0;
}
