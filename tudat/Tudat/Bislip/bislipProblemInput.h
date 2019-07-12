#ifndef BISLIPPROBLEMINPUT_H
#define BISLIPPROBLEMINPUT_H

#include <Tudat/Bislip/bislipParameters.h>


namespace bislip {

//! Wrapper class that contains the relevant hardware systems of a vehicle.
/*!
 *  Wrapper class that contains the relevant hardware systems of a vehicle. Not all member objects need to be set; nullptr
 *  member objects denote that the vehicle does not contain the associated hardware.
 */
class ProblemInput
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param playTimeString String containing a time-stamp of the entire batch (not defined; "" by default).
     */
    ProblemInput( const std::string playTimeString = "" ):
        playTimestring_( playTimeString ){ }

    //! Destructor
    ~ProblemInput( ){ }

    void setPlayTimePair( const std::pair< std::chrono::time_point< std::chrono::system_clock >, std::string > playTimePair )
    {
        playTimePair_ = playTimePair;
    }

   std::pair< std::chrono::time_point< std::chrono::system_clock >, std::string > getPlayTimePair( )
    {
        return playTimePair_;
    }

    void setProblemName( const std::string problemName )
    {
        problemName_ = problemName;
    }

    std::string getProblemName( )
    {
        return problemName_;
    }

    void setOutputPath( const std::string outputPath )
    {
        outputPath_ = outputPath;
    }

    std::string getOutputPath( )
    {
        return outputPath_;
    }

    void setOutputSubFolder( const std::string outputSubFolder )
    {
        outputSubFolder_ = outputSubFolder;
    }

    std::string getOutputSubFolder( )
    {
        return outputSubFolder_;
    }

    void setOutputSettings( const std::vector< double >  setOutputSettings )
    {
        setOutputSettings_ = setOutputSettings;
    }

    std::vector< double > getOutputSettings( )
    {
        return setOutputSettings_;
    }

    void setSimulationSettings( const std::vector< double > simulationSettings )
    {
        simulationSettings_ = simulationSettings;
    }

    std::vector< double > getSimulationSettings( )
    {
        return simulationSettings_;
    }

    void setVehicleName( const std::string vehicleName )
    {
        vehicleName_ = vehicleName;
    }

    std::string getVehicleName( )
    {
        return vehicleName_;
    }

    void setVehicleParameterList( const std::vector< std::string > vehicleParameterList )
    {
        vehicleParameterList_ = vehicleParameterList;
    }

    std::vector< std::string > getVehicleParameterList( )
    {
        return vehicleParameterList_;
    }

    void setVehicleParameters( const std::vector< double > vehicleParameters )
    {
        vehicleParameters_ = vehicleParameters;
    }

    std::vector< double > getVehicleParameters( )
    {
        return vehicleParameters_;
    }

    void setAerodynamicDatabaseFileList( const std::vector< std::string > aerodynamicDatabaseFileList )
    {
        aerodynamicDatabaseFileList_ = aerodynamicDatabaseFileList;
    }

    std::vector< std::string > getAerodynamicDatabaseFileList( )
    {
        return aerodynamicDatabaseFileList_;
    }

    void setInitialConditions( const std::vector< double > initialConditions )
    {
        initialConditions_ = initialConditions;
    }

    std::vector< double > getInitialConditions( )
    {
        return initialConditions_;
    }

    void setInitialState_Spherical( const Eigen::Vector6d initialState_Spherical )
    {
        initialState_Spherical_ = initialState_Spherical;
    }

    Eigen::Vector6d getInitialState_Spherical( )
    {
        return initialState_Spherical_;
    }

    void setAscentParameterList( const std::vector< std::string > ascentParameterList )
    {
        ascentParameterList_ = ascentParameterList;
    }

    std::vector< std::string > getAscentParameterList( )
    {
        return ascentParameterList_;
    }

    void setAscentParameterBounds( const std::vector< double > ascentParameterBounds )
    {
        ascentParameterBounds_ = ascentParameterBounds;
    }

    std::vector< double > getAscentParameterBounds( )
    {
        return ascentParameterBounds_;
    }

    void setAscentParameterBoundsMap( const std::map< bislip::Parameters::Interpolators, Eigen::MatrixXd > ascentParameterBoundsMap )
    {
        ascentParameterBoundsMap_ = ascentParameterBoundsMap;
    }

    const std::map< bislip::Parameters::Interpolators, Eigen::MatrixXd > getAscentParameterBoundsMap( )
    {
        return ascentParameterBoundsMap_;
    }

    void setDescentParameterList( const std::vector< std::string > descentParameterList )
    {
        descentParameterList_ = descentParameterList;
    }

    std::vector< std::string > getDescentParameterList( )
    {
        return descentParameterList_;
    }

    void setDescentParameterBounds( const std::vector< double > descentParameterBounds )
    {
        descentParameterBounds_ = descentParameterBounds;
    }

    std::vector< double > getDescentParameterBounds( )
    {
        return descentParameterBounds_;
    }

    void setDescentParameterBoundsMap( const std::map< bislip::Parameters::Interpolators, Eigen::MatrixXd > descentParameterBoundsMap )
    {
        descentParameterBoundsMap_ = descentParameterBoundsMap;
    }

    const std::map< bislip::Parameters::Interpolators, Eigen::MatrixXd > getDescentParameterBoundsMap( )
    {
        return descentParameterBoundsMap_;
    }

    void setDecisionVectorBounds( const std::vector< std::vector< double > > decisionVectorBounds )
    {
        decisionVectorBounds_ = decisionVectorBounds;
    }

    std::vector< std::vector< double > > getDecisionVectorBounds( )
    {
        return decisionVectorBounds_;
    }

    void setCentralBodies( const std::vector< std::string > centralBodies )
    {
        centralBodies_ = centralBodies;
    }

    std::vector< std::string > getCentralBodies( )
    {
        return centralBodies_;
    }

    void setBodiesToIntegrate( const std::vector< std::string > bodiesToIntegrate )
    {
        bodiesToIntegrate_ = bodiesToIntegrate;
    }

    std::vector< std::string > getBodiesToIntegrate( )
    {
        return bodiesToIntegrate_;
    }

    void setEarthRotationalEphemeris( const std::shared_ptr< tudat::ephemerides::RotationalEphemeris > earthRotationalEphemeris )
    {
        earthRotationalEphemeris_ = earthRotationalEphemeris;
    }

    std::shared_ptr< tudat::ephemerides::RotationalEphemeris > getEarthRotationalEphemeris( )
    {
        return earthRotationalEphemeris_;
    }

    void setAccelerationMap( const tudat::basic_astrodynamics::AccelerationMap accelerationMap )
    {
        accelerationMap_ = accelerationMap;
    }

    tudat::basic_astrodynamics::AccelerationMap getAccelerationMap( )
    {
        return accelerationMap_;
    }

    void setAscentTerminationSettings( const std::shared_ptr< tudat::propagators::PropagationTerminationSettings > ascentTerminationSettings )
    {
        ascentTerminationSettings_ = ascentTerminationSettings;
    }

    std::shared_ptr< tudat::propagators::PropagationTerminationSettings > getAscentTerminationSettings( )
    {
        return ascentTerminationSettings_;
    }

    void setDescentTerminationSettings( const std::shared_ptr< tudat::propagators::PropagationTerminationSettings > descentTerminationSettings )
    {
        descentTerminationSettings_ = descentTerminationSettings;
    }

    std::shared_ptr< tudat::propagators::PropagationTerminationSettings > getDescentTerminationSettings( )
    {
        return descentTerminationSettings_;
    }

    void setDependentVariablesToSave( const std::shared_ptr< tudat::propagators::DependentVariableSaveSettings > dependentVariablesToSave )
    {
        dependentVariablesToSave_ = dependentVariablesToSave;
    }

    std::shared_ptr< tudat::propagators::DependentVariableSaveSettings > getDependentVariablesToSave( )
    {
        return dependentVariablesToSave_;
    }

    void setConstraints( const std::vector< double > constraintsValues )
    {
        constraintsValues_ = constraintsValues;
    }

    std::vector< double > getConstraints( )
    {
        return constraintsValues_;
    }

    void setHardConstraints( const std::vector< double > hardConstraintsValues )
    {
        hardConstraintsValues_ = hardConstraintsValues;
    }

    std::vector< double > getHardConstraints( )
    {
        return hardConstraintsValues_;
    }


    void setPopulation( const  std::map < std::string, Eigen::VectorXd > population )
    {
        population_ = population;
    }

    std::map < std::string, Eigen::VectorXd > getPopulation( )
    {
        return population_;
    }

    void setFitness( const  std::map < std::string, Eigen::VectorXd > fitness )
    {
        fitness_ = fitness;
    }

    std::map < std::string, Eigen::VectorXd > getFitness( )
    {
        return fitness_;
    }

    void setPrintedPopulation( const  std::map < std::string, Eigen::VectorXd > printedPopulation )
    {
        printedPopulation_ = printedPopulation;
    }

    std::map < std::string, Eigen::VectorXd > getPrintedPopulation( )
    {
        return printedPopulation_;
    }

    void setPrintedFitness( const  std::map < std::string, Eigen::VectorXd > printedFitness )
    {
        printedFitness_ = printedFitness;
    }

    std::map < std::string, Eigen::VectorXd > getPrintedFitness( )
    {
        return printedFitness_;
    }
    void setIndividualNumber( const int individualNumber )
    {
        individualNumber_ = individualNumber;
    }

    int getIndividualNumber( )
    {
        return individualNumber_;
    }
    void setObjectiveFunctionCase( const char objectiveFunctionCase )
    {
        objectiveFunctionCase_ = objectiveFunctionCase;
    }

    char getObjectiveFunctionCase( )
    {
        return objectiveFunctionCase_;
    }


/*
    void setAtmosphericModel_US76( const tudat::aerodynamics::TabulatedAtmosphere atmosphericModel_US76 )
    {
        atmosphericModel_US76_ = atmosphericModel_US76;
    }

    tudat::aerodynamics::TabulatedAtmosphere getAtmosphericModel_US76( )
    {
        return atmosphericModel_US76_;
    }

    void setAtmosphericModel_NRLMSISE00( const tudat::aerodynamics::NRLMSISE00Atmosphere atmosphericModel_NRLMSISE00 )
    {
        atmosphericModel_NRLMSISE00_ = atmosphericModel_NRLMSISE00;
    }

    tudat::aerodynamics::NRLMSISE00Atmosphere getAtmosphericModel_NRLMSISE00( )
    {
        return atmosphericModel_NRLMSISE00_;
    }
*/
    void setEvolutionEvaluationFlag( const bool evolutionEvaluationFlag )
    {
        evolutionEvaluationFlag_ = evolutionEvaluationFlag;
    }

    bool getEvolutionEvaluationFlag( )
    {
        return evolutionEvaluationFlag_;
    }
    void setRerunFileNameSuffix( const std::string rerunFileNameSuffix )
    {
        rerunFileNameSuffix_ = rerunFileNameSuffix;
    }

    std::string getRerunFileNameSuffix( )
    {
        return rerunFileNameSuffix_;
    }





private:
    //tudat::aerodynamics::NRLMSISE00Atmosphere atmosphericModel_NRLMSISE00_;
    //tudat::aerodynamics::TabulatedAtmosphere atmosphericModel_US76_;

    std::string rerunFileNameSuffix_;
    bool evolutionEvaluationFlag_;
    char objectiveFunctionCase_;
    int individualNumber_;
    std::map < std::string, Eigen::VectorXd > population_;
    std::map < std::string, Eigen::VectorXd > fitness_;
    std::map < std::string, Eigen::VectorXd > printedPopulation_;
    std::map < std::string, Eigen::VectorXd > printedFitness_;
    std::pair< std::chrono::time_point< std::chrono::system_clock >, std::string > playTimePair_;
    std::string playTimestring_;
    std::string problemName_;
    std::string outputPath_;
    std::string outputSubFolder_;
    std::vector< double > hardConstraintsValues_;
    std::vector< double > setOutputSettings_;
    std::vector< double > simulationSettings_;
    std::string vehicleName_;
    std::vector< std::string > vehicleParameterList_;
    std::vector< double > vehicleParameters_;
    std::vector< std::string > aerodynamicDatabaseFileList_;
    std::vector< double > initialConditions_;
    Eigen::Vector6d initialState_Spherical_;
    std::vector< std::string > ascentParameterList_;
    std::vector< double > ascentParameterBounds_;
    std::map< bislip::Parameters::Interpolators, Eigen::MatrixXd > ascentParameterBoundsMap_;
    std::vector< std::string > descentParameterList_;
    std::vector< double > descentParameterBounds_;
    std::map< bislip::Parameters::Interpolators, Eigen::MatrixXd > descentParameterBoundsMap_;
    std::vector< std::vector< double > > decisionVectorBounds_;
    std::vector< std::string > centralBodies_;
    std::vector< std::string > bodiesToIntegrate_;
    std::shared_ptr< tudat::ephemerides::RotationalEphemeris > earthRotationalEphemeris_;
    tudat::basic_astrodynamics::AccelerationMap accelerationMap_;
    std::shared_ptr< tudat::propagators::PropagationTerminationSettings > ascentTerminationSettings_;
    std::shared_ptr< tudat::propagators::PropagationTerminationSettings > descentTerminationSettings_;
    std::shared_ptr< tudat::propagators::DependentVariableSaveSettings > dependentVariablesToSave_;
    std::vector< double > constraintsValues_;

};

} // namespace bislip

#endif // BISLIPPROBLEMINPUT_H
