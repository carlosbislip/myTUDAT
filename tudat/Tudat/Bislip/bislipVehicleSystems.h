#ifndef BISLIPVEHICLESYSTEMS_H
#define BISLIPVEHICLESYSTEMS_H

#include <Tudat/Bislip/bislipParameters.h>


namespace bislip {

//! Wrapper class that contains the relevant hardware systems of a vehicle.
/*!
 *  Wrapper class that contains the relevant hardware systems of a vehicle. Not all member objects need to be set; nullptr
 *  member objects denote that the vehicle does not contain the associated hardware.
 */
class BislipVehicleSystems
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param dryMass Total dry mass of the vehicle (not defined; NaN by default).
     */
    BislipVehicleSystems( const double landingMass = 0 ):
        landingMass_( landingMass ){ }

    //! Destructor
    ~BislipVehicleSystems( ){ }


    void setInitialMass( const double initialMass )
    {
        initialMass_ = initialMass;
    }

    double getLandingMass( )
    {
        return landingMass_;
    }

    //! Function to retrieve the initial mass of the vehicle
    /*!
     * Function to retrieve the  total dry mass of the vehicle
     * \return initial mass of the vehicle
     */
    double getInitialMass( )
    {
        return initialMass_;
    }

    void setDebugInfo( const int debugInfo )
    {
        debugInfo_ = debugInfo;
    }

    int getDebugInfo( )
    {
        return debugInfo_;
    }
    void setAverageEarthRadius( const double averageEarthRadius )
    {
        averageEarthRadius_ = averageEarthRadius;
    }

    double getAverageEarthRadius( )
    {
        return averageEarthRadius_;
    }

    void setReferenceArea( const double referenceArea )
    {
        referenceArea_ = referenceArea;
    }

    double getReferenceArea( )
    {
        return referenceArea_;
    }

    void setMaxThrust( const double maxThrust )
    {
        maxThrust_ = maxThrust;
    }

    double getMaxThrust( )
    {
        return maxThrust_;
    }

    void setSpecificImpulse( const double specificImpulse )
    {
        specificImpulse_ = specificImpulse;
    }

    double getSpecificImpulse( )
    {
        return specificImpulse_;
    }

    void setNoseRadius( const double noseRadius )
    {
        noseRadius_ = noseRadius;
    }

    double getNoseRadius( )
    {
        return noseRadius_;
    }

    void setLeadingEdgeRadius( const double leadingEdgeRadius )
    {
        leadingEdgeRadius_ = leadingEdgeRadius;
    }

    double getLeadingEdgeRadius( )
    {
        return leadingEdgeRadius_;
    }
    void setWorkingRadius( const double workingRadius )
    {
        workingRadius_ = workingRadius;
    }

    double getWorkingRadius( )
    {
        return workingRadius_;
    }

    void setNoseOrLeadingEdge( const bool noseOrLeadingEdge )
    {
        noseOrLeadingEdge_ = noseOrLeadingEdge;
    }

    bool getNoseOrLeadingEdge( )
    {
        return noseOrLeadingEdge_;
    }



    void setWingSweepAngle( const double lambda )
    {
        lambda_ = lambda;
    }

    double getWingSweepAngle( )
    {
        return lambda_;
    }

    void setLocalBodyAngle( const double phi )
    {
        phi_ = phi;
    }

    double getLocalBodyAngle( )
    {
        return phi_;
    }

    void setTransitionDistance( const double x_T )
    {
        x_T_ = x_T;
    }

    double getTransitionDistance( )
    {
        return x_T_;
    }

    void setE_max( const double E_max )
    {
        E_max_ = E_max;
    }

    double getE_max( )
    {
        return E_max_;
    }
    void setMinimumDynamicPressureforControlSurface( const double minimumDynamicPressureforControlSurface )
    {
        minimumDynamicPressureforControlSurface_ = minimumDynamicPressureforControlSurface;
    }

    double getMinimumDynamicPressureforControlSurface( )
    {
        return minimumDynamicPressureforControlSurface_;
    }

    void setBodyFlapDeflectionLimits( const std::pair < double, double > bodyFlapBounds )
    {
        bodyFlapBounds_ = bodyFlapBounds;
    }

    std::pair < double, double > getBodyFlapDeflectionLimits( )
    {
        return bodyFlapBounds_;
    }
    void setElevonDeflectionLimits( const std::pair < double, double > elevonBounds )
    {
        elevonBounds_ = elevonBounds;
    }

    std::pair < double, double > getElevonDeflectionLimits( )
    {
        return elevonBounds_;
    }






    /*
    void setAngleOfAttackInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack )
    {
        interpolatorAngleOfAttack_ = interpolatorAngleOfAttack;
    }
    void setBankAngleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngle )
    {
        interpolatorBankAngle_ = interpolatorBankAngle;
    }
    void setThrustElevationAngleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngle )
    {
        interpolatorThrustElevationAngle_ = interpolatorThrustElevationAngle;
    }
    void setThrustAzimuthAngleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle )
    {
        interpolatorThrustAzimuthAngle_ = interpolatorThrustAzimuthAngle;
    }
    void setThrottleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting )
    {
        interpolatorThrottleSetting_ = interpolatorThrottleSetting;
    }

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getAngleOfAttackInterpolator( )
    {
        return interpolatorAngleOfAttack_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getBankAngleInterpolator( )
    {
        return interpolatorBankAngle_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrustElevationAngleInterpolator( )
    {
        return interpolatorThrustElevationAngle_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrustAzimuthAngleInterpolator( )
    {
        return interpolatorThrustAzimuthAngle_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrottleInterpolator( )
    {
        return interpolatorThrottleSetting_;
    }
*/
    void setTargetLat( const double target_lat )
    {
        target_lat_ = target_lat;
    }

    void setTargetLon( const double target_lon )
    {
        target_lon_ = target_lon;
    }

    void setInitialLat( const double initial_lat )
    {
        initial_lat_ = initial_lat;
    }
    void setInitialLon( const double initial_lon )
    {
        initial_lon_ = initial_lon;
    }
    void setInitialDistanceToTarget( const double initial_d_to_target )
    {
        initial_d_to_target_ = initial_d_to_target;
    }
    void setFinalDistanceToTarget( const double final_d_to_target )
    {
        final_d_to_target_ = final_d_to_target;
    }
    double getTargetLat( )
    {
        return target_lat_;
    }
    double getTargetLon( )
    {
        return target_lon_;
    }
    double getInitialLat( )
    {
        return initial_lat_;
    }
    double getInitialLon( )
    {
        return initial_lon_;
    }
    double getInitialDistanceToTarget( )
    {
        return initial_d_to_target_;
    }
    double getFinalDistanceToTarget( )
    {
        return final_d_to_target_;
    }
    void setParameterBounds( const std::map< bislip::Parameters::Optimization, std::pair < double, double > > &Bounds )
    {
        Bounds_ = Bounds;
    }
    std::pair < double, double > getParameterBounds( const bislip::Parameters::Optimization &parameter )
    {
        return Bounds_[ parameter ];
    }
    void setParameterInterpolator( const std::map< bislip::Parameters::Optimization, std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > > &Interpolators )
    {
        Interpolators_ = Interpolators;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getParameterInterpolator( const bislip::Parameters::Optimization &parameter )
    {
        return Interpolators_[ parameter ];
    }
    void setInitialCoordinates( const std::pair < double, double > &initialCoordinates )
    {
        initialCoordinates_ = initialCoordinates;
    }
    std::pair < double, double > getInitialCoordinates( )
    {
        return initialCoordinates_;
    }
    void setTargetCoordinates( const std::pair < double, double > &targetCoordinates )
    {
        targetCoordinates_ = targetCoordinates;
    }
    std::pair < double, double > getTargetCoordinates( )
    {
        return targetCoordinates_;
    }
    void setStartingEpoch( const double startingEpoch )
    {
        startingEpoch_ = startingEpoch;
    }
    double getStartingEpoch( )
    {
        return startingEpoch_;
    }
    void setPropagationStepSize( const double propagationStepSize )
    {
        propagationStepSize_ = propagationStepSize;
    }
    double getPropagationStepSize( )
    {
        return propagationStepSize_;
    }
    void setGuidanceStepSize( const double guidanceStepSize )
    {
        guidanceStepSize_ = guidanceStepSize;
    }
    double getGuidanceStepSize( )
    {
        return guidanceStepSize_;
    }

    void setSamplingRatio( const double samplingRatio )
    {
        samplingRatio_ = samplingRatio;
    }
    double getSamplingRatio( )
    {
        return samplingRatio_;
    }
    void setReferenceValues( const Eigen::Vector3d referenceValues )
    {
        referenceValues_ = referenceValues;
    }
    Eigen::Vector3d getReferenceValues( )
    {
        return referenceValues_;
    }
    void setMomentReferenceCenter( const Eigen::Vector3d momentReferenceCenter )
    {
        momentReferenceCenter_ = momentReferenceCenter;
    }
    Eigen::Vector3d getMomentReferenceCenter( )
    {
        return momentReferenceCenter_;
    }
    void setMassReferenceCenter( const Eigen::Vector3d massReferenceCenter )
    {
        massReferenceCenter_ = massReferenceCenter;
    }
    Eigen::Vector3d getMassReferenceCenter( )
    {
        return massReferenceCenter_;
    }
    void setThrustReferenceCenter( const Eigen::Vector3d thrustReferenceCenter )
    {
        thrustReferenceCenter_ = thrustReferenceCenter;
    }
    Eigen::Vector3d getThrustReferenceCenter( )
    {
        return thrustReferenceCenter_;
    }

    void setAlphaMachEnvelopeUBInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAlphaMachEnvelopeUB )
    {
        interpolatorAlphaMachEnvelopeUB_ = interpolatorAlphaMachEnvelopeUB;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getAlphaMachEnvelopeUBInterpolator( )
    {
        return interpolatorAlphaMachEnvelopeUB_;
    }
    void setAlphaMachEnvelopeLBInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAlphaMachEnvelopeLB )
    {
        interpolatorAlphaMachEnvelopeLB_ = interpolatorAlphaMachEnvelopeLB;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getAlphaMachEnvelopeLBInterpolator( )
    {
        return interpolatorAlphaMachEnvelopeLB_;
    }
    void setHeadingErrorDeadBandInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorHeadingErrorDeadband )
    {
        interpolatorHeadingErrorDeadband_ = interpolatorHeadingErrorDeadband;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getHeadingErrorDeadBandInterpolator( )
    {
        return interpolatorHeadingErrorDeadband_;
    }


    void setKourouAngleOfAttackInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorKourouAngleOfAttack )
    {
        interpolatorKourouAngleOfAttack_ = interpolatorKourouAngleOfAttack;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getKourouAngleOfAttackInterpolator( )
    {
        return interpolatorKourouAngleOfAttack_;
    }


    void setKourouBankAngleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorKourouBankAngle )
    {
        interpolatorKourouBankAngle_ = interpolatorKourouBankAngle;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getKourouBankAngleInterpolator( )
    {
        return interpolatorKourouBankAngle_;
    }
    void setValidationFlag( const bool validationFlag )
    {
        validationFlag_ = validationFlag;
    }
    bool getValidationFlag( )
    {
        return validationFlag_;
    }


    void setCurrentTrajectoryPhase( const std::string trajectoryPhase )
    {
        trajectoryPhase_ = trajectoryPhase;
    }
    std::string getCurrentTrajectoryPhase( )
    {
        return trajectoryPhase_;
    }
    void setCurrentAngleOfAttack( const double currentAngleOfAttack )
    {
        currentAngleOfAttack_ = currentAngleOfAttack;
    }
    double getCurrentAngleOfAttack( )
    {
        return currentAngleOfAttack_;
    }
    void setCurrentBodyFlapAngle( const double currentBodyFlapAngle )
    {
        currentBodyFlapAngle_ = currentBodyFlapAngle;
    }
    double getCurrentBodyFlapAngle( )
    {
        return currentBodyFlapAngle_;
    }
    void setCurrentElevonAngle( const double currentElevonAngle )
    {
        currentElevonAngle_ = currentElevonAngle;
    }
    double getCurrentElevonAngle( )
    {
        return currentElevonAngle_;
    }

    void setEvaluatedBankAngle( const double evaluatedBankAngle )
    {
        evaluatedBankAngle_ = evaluatedBankAngle;
    }
    double getEvaluatedBankAngle( )
    {
        return evaluatedBankAngle_;
    }

    void setCurrentBankAngle( const double currentBankAngle )
    {
        currentBankAngle_ = currentBankAngle;
    }
    double getCurrentBankAngle( )
    {
        return currentBankAngle_;
    }
    void setTempBankAngle( const double tempBankAngle )
    {
        tempBankAngle_ = tempBankAngle;
    }
    double getTempBankAngle( )
    {
        return tempBankAngle_;
    }
    void setCurrentBankAngleLimit( const double currentBankAngleLimit )
    {
        currentBankAngleLimit_ = currentBankAngleLimit;
    }
    double getCurrentBankAngleLimit( )
    {
        return currentBankAngle_;
    }
    void setBankAngleReversalTrigger( const int bankAngleReversalTrigger )
    {
        bankAngleReversalTrigger_ = bankAngleReversalTrigger;
    }
    int getBankAngleReversalTrigger( )
    {
        return bankAngleReversalTrigger_;
    }
    void setBankAngleReversalTimepoint( const double bankAngleReversalTimepoint )
    {
        bankAngleReversalTimepoint_ = bankAngleReversalTimepoint;
    }
    double getBankAngleReversalTimepoint( )
    {
        return bankAngleReversalTimepoint_;
    }


    void setRelativeChestForwardAngle( const double relativeChestForwardAngle )
    {
        relativeChestForwardAngle_ = relativeChestForwardAngle;
    }
    double getRelativeChestForwardAngle( )
    {
        return relativeChestForwardAngle_;
    }
    void setVertebralColumnInclinationAngle( const double vertebralColumnInclinationAngle )
    {
        vertebralColumnInclinationAngle_ = vertebralColumnInclinationAngle;
    }
    double getVertebralColumnInclinationAngle( )
    {
        return vertebralColumnInclinationAngle_;
    }

    void setBodyFrameToPassengerFrameTransformationMatrix( const Eigen::Matrix3d bodyFrameToPassengerFrameTransformationMatrix )
    {
        bodyFrameToPassengerFrameTransformationMatrix_ = bodyFrameToPassengerFrameTransformationMatrix;
    }
    Eigen::Matrix3d getBodyFrameToPassengerFrameTransformationMatrix( )
    {
        return bodyFrameToPassengerFrameTransformationMatrix_;
    }
    void setCurrentBodyFixedThrustDirection( const Eigen::Vector3d currentBodyFixedThrustDirection )
    {
        currentBodyFixedThrustDirection_ = currentBodyFixedThrustDirection;
    }
    Eigen::Vector3d getCurrentBodyFixedThrustDirection( )
    {
        return currentBodyFixedThrustDirection_;
    }
    void setCurrentEngineStatus( const bool currentEngineStatus )
    {
        currentEngineStatus_ = currentEngineStatus;
    }
    bool getCurrentEngineStatus( )
    {
        return currentEngineStatus_;
    }
    void setCurrentThrustMagnitude( const double currentThrustMagnitude )
    {
        currentThrustMagnitude_ = currentThrustMagnitude;
    }
    double getCurrentThrustMagnitude( )
    {
        return currentThrustMagnitude_;
    }
    void setThrustMagnitudeOutput( const double thrustMagnitudeOutput )
    {
        thrustMagnitudeOutput_ = thrustMagnitudeOutput;
    }
    double getThrustMagnitudeOutput( )
    {
        return thrustMagnitudeOutput_;
    }
    void setBodyFixedThrustVectorOutput( const Eigen::Vector3d bodyFixedThrustVectorOutput )
    {
        bodyFixedThrustVectorOutput_ = bodyFixedThrustVectorOutput;
    }
    Eigen::Vector3d getBodyFixedThrustVectorOutput( )
    {
        return bodyFixedThrustVectorOutput_;
    }
    void setCurrentThrustElevationAngle( const double currentThrustElevationAngle )
    {
        currentThrustElevationAngle_ = currentThrustElevationAngle;
    }
    double getCurrentThrustElevationAngle( )
    {
        return currentThrustElevationAngle_;
    }
    void setCurrentThrustAzimuthAngle( const double currentThrustAzimuthAngle )
    {
        currentThrustAzimuthAngle_ = currentThrustAzimuthAngle;
    }
    double getCurrentThrustAzimuthAngle( )
    {
        return currentThrustAzimuthAngle_;
    }
    void setCurrentThrottleSetting( const double currentThrottleSetting )
    {
        currentThrottleSetting_ = currentThrottleSetting;
    }
    double getCurrentThrottleSetting( )
    {
        return currentThrottleSetting_;
    }
    void setReversalConditional( const double reversalConditional )
    {
        reversalConditional_ = reversalConditional;
    }
    double getReversalConditional( )
    {
        return reversalConditional_;
    }

    void setWallTemperature( const double wallTemperature )
    {
        wallTemperature_ = wallTemperature;
    }
    double getWallTemperature( )
    {
        return wallTemperature_;
    }

    void setCurrentHeatFluxStagnation( const double currentHeatFluxStagnation )
    {
        currentHeatFluxStagnation_ = currentHeatFluxStagnation;
    }
    double setCurrentHeatFluxStagnation( )
    {
        return currentHeatFluxStagnation_;
    }

    void setCurrentHeatFluxTauber( const double currentHeatFluxTauber )
    {
        currentHeatFluxTauber_ = currentHeatFluxTauber;
    }
    double getCurrentHeatFluxTauber( )
    {
        return currentHeatFluxTauber_;
    }
    void setTauberHeatFluxStagnation( const double tauberHeatFluxStagnation )
    {
        tauberHeatFluxStagnation_ = tauberHeatFluxStagnation;
    }
    double getTauberHeatFluxStagnation( )
    {
        return tauberHeatFluxStagnation_;
    }

    void setTauberHeatFluxFlatPlate( const double tauberHeatFluxFlatPlate )
    {
        tauberHeatFluxFlatPlate_ = tauberHeatFluxFlatPlate;
    }
    double getTauberHeatFluxFlatPlate( )
    {
        return tauberHeatFluxFlatPlate_;
    }

    void setCurrentHeatFluxChapman( const double currentHeatFluxChapman )
    {
        currentHeatFluxChapman_ = currentHeatFluxChapman;
    }
    double getCurrentHeatFluxChapman( )
    {
        return currentHeatFluxChapman_;
    }
    void setChapmanWallTemp( const double chapmanWallTemp )
    {
        chapmanWallTemp_ = chapmanWallTemp;
    }
    double getChapmanWallTemp( )
    {
        return chapmanWallTemp_;
    }
    void setTauberWallTempStagnation( const double tauberWallTempStagnation )
    {
        tauberWallTempStagnation_ = tauberWallTempStagnation;
    }
    double getTauberWallTempStagnation( )
    {
        return tauberWallTempStagnation_;
    }
    void setTauberWallTempFlatPlate( const double tauberWallTempFlatPlate )
    {
        tauberWallTempFlatPlate_ = tauberWallTempFlatPlate;
    }
    double getTauberWallTempFlatPlate( )
    {
        return tauberWallTempFlatPlate_;
    }

    void setPreviousCoordinates( const Eigen::Vector2d previousCoordinates )
    {
        previousCoordinates_ = previousCoordinates;
    }
    Eigen::Vector2d getPreviousCoordinates( )
    {
        return previousCoordinates_;
    }

    void setPreviousCartesianCoordinates( const Eigen::Vector3d previousCartesianCoordinates )
    {
        previousCartesianCoordinates_ = previousCartesianCoordinates;
    }
    Eigen::Vector3d getPreviousCartesianCoordinates( )
    {
        return previousCartesianCoordinates_;
    }
    void setCumulativeDistanceTravelled( const double cumulativeDistanceTravelled )
    {
        cumulativeDistanceTravelled_ = cumulativeDistanceTravelled;
    }
    double getCumulativeDistanceTravelled( )
    {
        return cumulativeDistanceTravelled_;
    }
    void setCumulativeAngularDistanceTravelled( const double groundtrackDistance )
    {
        groundtrackDistance_ = groundtrackDistance;
    }
    double getCumulativeAngularDistanceTravelled( )
    {
        return groundtrackDistance_;
    }
    void setCurrentFlightPathAngleRate( const double currentFlightPathAngleRate )
    {
        currentFlightPathAngleRate_ = currentFlightPathAngleRate;
    }
    double getCurrentFlightPathAngleRate( )
    {
        return currentFlightPathAngleRate_;
    }
    void setSkipSuppressionTimingTrigger( const double skipSuppressionTimingTrigger )
    {
        skipSuppressionTimingTrigger_ = skipSuppressionTimingTrigger;
    }
    double getSkipSuppressionTimingTrigger( )
    {
        return skipSuppressionTimingTrigger_;
    }

    void setAscentTerminationDistanceRatio( const double ascentTerminationDistanceRatio )
    {
        ascentTerminationDistanceRatio_ = ascentTerminationDistanceRatio;
    }
    double getAscentTerminationDistanceRatio( )
    {
        return ascentTerminationDistanceRatio_;
    }
    void setInitialAltitude( const double initialAltitude )
    {
        initialAltitude_ = initialAltitude;
    }
    double getInitialAltitude( )
    {
        return initialAltitude_;
    }
    void setInitialAirspeed( const double initialAirspeed )
    {
        initialAirspeed_ = initialAirspeed;
    }
    double getInitialAirspeed( )
    {
        return initialAirspeed_;
    }


private:

    double initialAltitude_;
    double initialAirspeed_;
    double ascentTerminationDistanceRatio_;
    double skipSuppressionTimingTrigger_;
    std::pair< double, double > bodyFlapBounds_;
    std::pair< double, double > elevonBounds_;

    Eigen::Vector3d previousCartesianCoordinates_;
    Eigen::Vector2d previousCoordinates_;
    double cumulativeDistanceTravelled_;
    double currentFlightPathAngleRate_;
    double groundtrackDistance_;
    double currentHeatFluxChapman_;
    double chapmanWallTemp_;
    double tauberWallTempStagnation_;
    double tauberWallTempFlatPlate_;
    double averageEarthRadius_;

    int debugInfo_;
    double currentAngleOfAttack_;
    double currentBodyFlapAngle_;
    double currentElevonAngle_;
    double evaluatedBankAngle_;
    double minimumDynamicPressureforControlSurface_;
    double currentBankAngle_;
    double reversalConditional_;
    double tempBankAngle_;
    double currentBankAngleLimit_;
    double currentThrustMagnitude_;
    double currentThrustAzimuthAngle_;
    double currentThrustElevationAngle_;
    double currentThrottleSetting_;
    Eigen::Vector3d currentBodyFixedThrustDirection_;
    bool currentEngineStatus_;
    double thrustMagnitudeOutput_;
    Eigen::Vector3d bodyFixedThrustVectorOutput_;
    int bankAngleReversalTrigger_;
    double bankAngleReversalTimepoint_;



    double relativeChestForwardAngle_;
    double vertebralColumnInclinationAngle_;
    Eigen::Matrix3d bodyFrameToPassengerFrameTransformationMatrix_;

    std::string trajectoryPhase_;

    std::pair < double, double > initialCoordinates_;
    std::pair < double, double > targetCoordinates_;

    double workingRadius_;
    double noseRadius_;
    double leadingEdgeRadius_;
    bool noseOrLeadingEdge_;
    double lambda_;
    double phi_;
    double x_T_;
    double E_max_;
    double maxThrust_;
    double specificImpulse_;
    double startingEpoch_;
    double propagationStepSize_;
    double guidanceStepSize_;
    double samplingRatio_;
    double referenceArea_;
    // double initialMass_;

    double wallTemperature_;
    double currentHeatFluxStagnation_;
    double currentHeatFluxTauber_;
    double tauberHeatFluxStagnation_;
    double tauberHeatFluxFlatPlate_;

    //    std::map< std::string, std::pair < double, double > > Bounds_;
    // std::map< std::string, std::pair < double, double > > Bounds_Ascent_;
    //std::map< std::string, std::pair < double, double > > Bounds_Descent_;
    std::map< bislip::Parameters::Optimization, std::pair < double, double > > Bounds_;
    std::map< bislip::Parameters::Optimization, std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > > Interpolators_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAlphaMachEnvelopeUB_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAlphaMachEnvelopeLB_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorHeadingErrorDeadband_;


std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorKourouAngleOfAttack_;
std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorKourouBankAngle_;
bool validationFlag_;
    /*
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngle_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngle_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting_;


    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngleAscent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngleAscent_;

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngleDescent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngleDescent_;
   */

    //! Passing target coordinates
    double target_lat_;
    double target_lon_;

    //! Passing initial coordinates
    double initial_lat_;
    double initial_lon_;

    //! Passing initial distance to target
    double initial_d_to_target_;

    //! Passing final distance to target
    double final_d_to_target_;

    //! Initial mass of the vehicle
    double initialMass_;

    //! Landing mass of the vehicle
    double landingMass_;


    //double b_ref_;
    //double c_ref_;
    //double del_x_;
    //double del_x_T_;
    //double del_z_T_;
    Eigen::Vector3d referenceValues_;
    Eigen::Vector3d momentReferenceCenter_;
    Eigen::Vector3d thrustReferenceCenter_;
    Eigen::Vector3d massReferenceCenter_;

};

} // namespace bislip

//pagmo::algorithm getPagmoAlgorithm ( const int index );

#endif // BISLIPVEHICLESYSTEMS_H
