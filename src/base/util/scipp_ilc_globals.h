namespace scipp_ilc {

    //Various sidloi3-IR_realign geometric constants.
    //All measurements in units of millimeters.
    static const float _LumiCal_zmin = 1557.0;
    static const float _LumiCal_thickness = 138.5;

    static const float _BeamCal_zmin = 3265;
    static const float _BeamCal_thickness = 175.0;
    static const float _BeamCal_outer_radius = 140.0;
    static const float _BeamCal_outgoing_pipe_radius = 20.5;
    static const float _BeamCal_incoming_pipe_radius = 15.5;

    static const float _crossing_angle = 0.014; // radians
    static const float _transform = _crossing_angle / 2.0;



    //I am removing everything which hits at the very edge of the beamcal
    //(two 3.5 mm pixels from the 140 mm edge) in order to deal with the
    //poor statistical data at the outer boundries
    static const float _radius_cut = _BeamCal_outer_radius - 7.0; //mm.
}
