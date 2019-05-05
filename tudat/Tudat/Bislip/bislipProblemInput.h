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
     * \param dryMass Total dry mass of the vehicle (not defined; NaN by default).
     */
    ProblemInput( const std::string playTime = "" ):
        playTime_( playTime ){ }

    //! Destructor
    ~ProblemInput( ){ }

    void setPlayTime( const std::string playTime )
    {
        playTime_ = playTime;
    }

   std::string getPlayTime( )
    {
        return playTime_;
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

    void setAscentParameterBoundsMap( const std::map< bislip::Parameters::Optimization, std::pair < double, double > > ascentParameterBoundsMap )
    {
        ascentParameterBoundsMap_ = ascentParameterBoundsMap;
    }

    const std::map< bislip::Parameters::Optimization, std::pair < double, double > > getAscentParameterBoundsMap( )
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

    void setDescentParameterBoundsMap( const std::map< bislip::Parameters::Optimization, std::pair < double, double > > descentParameterBoundsMap )
    {
        descentParameterBoundsMap_ = descentParameterBoundsMap;
    }

    const std::map< bislip::Parameters::Optimization, std::pair < double, double > > getDescentParameterBoundsMap( )
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


private:

    std::string playTime_;
    std::string problemName_;
    std::string outputPath_;
    std::string outputSubFolder_;
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
    std::map< bislip::Parameters::Optimization, std::pair < double, double > > ascentParameterBoundsMap_;
    std::vector< std::string > descentParameterList_;
    std::vector< double > descentParameterBounds_;
    std::map< bislip::Parameters::Optimization, std::pair < double, double > > descentParameterBoundsMap_;
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
