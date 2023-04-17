template_text ="""<?xml version="1.0" encoding="UTF-8"?>
<Plant name="Heliantus_Pagès_2013_new" filetype="parameters">
<organ type="seed" name="true" subType="0">
    <parameter name="seedPos.x" value="0"/>
    <parameter name="seedPos.y" value="0"/>
    <parameter name="seedPos.z" value="-3"/>
    <parameter name="plantingdepth" value="3"/>
    <parameter name="firstB" value="0" dev="0"/>
    <parameter name="delayB" value="0" dev="0"/>
    <parameter name="maxB" value="0" dev="0"/>
    <parameter name="maxTi" value="0" dev="0"/>
    <parameter name="nC" value="0"/>
    <parameter name="nz" value="0"/>
    <parameter name="firstSB" value="1000" dev="0"/>
    <parameter name="delaySB" value="1000" dev="0"/>
    <parameter name="delayRC" value="1000" dev="0"/>
    <parameter name="simulationTime" value="20" dev="0"/>
</organ>

<organ type="root" name="taproot" subType="1">
    <parameter name="lb" value="1" dev="0"/>
    <!--Basal zone [cm]-->
    <parameter name="la" value="1" dev="0"/>
    <!--Apical zone [cm];-->
    <parameter name="ln" value="5" dev="0" functiontype="0"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="lmax" value="102" dev="0"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="nob" value="11.9" dev="0"/>
    <parameter name="r" value="2" dev="0.45"/>
    <parameter name="a" value="0.2" />
    <parameter name="RotBeta" value="1"/>
    <parameter name="BetaDev" value="0.2"/>
    <parameter name="InitBeta" value="0"/>
    <parameter name="tropismT" value="4"/>
    <parameter name="tropismN" value="18"/>
    <parameter name="tropismS" value="0.01"/>
    <parameter name="dx" value="0.25"/>
    <parameter name="theta" value="0" dev="0"/>
    <parameter name="rlt" value="1000000000" dev="0"/>
    <parameter name="gf" value="1"/>
    {root}
</organ>

<organ type="root" name="lateral1" subType="2">
    <parameter name="lb" value="8" dev="0"/>
    <!--Basal zone [cm]-->
    <parameter name="la" value="1" dev="0"/>
    <!--Apical zone [cm];-->
    <parameter name="ln" value="10" dev="0" functiontype="1"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="lmax" value="23" dev="0"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="r" value="1"/>
    <parameter name="a" value="0.3" />
    <parameter name="RotBeta" value="1"/>
    <parameter name="BetaDev" value="0"/>
    <parameter name="InitBeta" value="0"/>
    <parameter name="tropismT" value="1"/>
    <parameter name="tropismN" value="2"/>
    <parameter name="tropismS" value="0.1"/>
    <parameter name="dx" value="0.25"/>
    <parameter name="theta" value="0" dev="0"/>
    <parameter name="rlt" value="1000000000" dev="0"/>
    <parameter name="gf" value="1"/>
</organ>
<organ type="root" name="lateral2" subType="3">
    <parameter name="lb" value="8" dev="0"/>
    <!--Basal zone [cm]-->
    <parameter name="la" value="1" dev="0"/>
    <!--Apical zone [cm];-->
    <parameter name="ln" value="10" dev="0" functiontype="1"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="lmax" value="10" dev="0"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="r" value="1"/>
    <parameter name="a" value="0.3" />
    <parameter name="RotBeta" value="1"/>
    <parameter name="BetaDev" value="0"/>
    <parameter name="InitBeta" value="0"/>
    <parameter name="tropismT" value="1"/>
    <parameter name="tropismN" value="2"/>
    <parameter name="tropismS" value="1"/>
    <parameter name="dx" value="0.1"/>
    <parameter name="theta" value="0" dev="0"/>
    <parameter name="rlt" value="1000000000" dev="0"/>
    <parameter name="gf" value="1"/>
</organ>


<organ type="stem" name="mainstem" subType="1">
    <parameter name="lb" value="1" dev="0"/>
    <!--Basal zone [cm]-->
    <parameter name="la" value="1" dev="0"/>
    <!--Apical zone [cm];-->
    <parameter name="ln" value="10" dev="0" functiontype="0"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="lmax" value="102" dev="0"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="nob" value="11.9" dev="0"/>
    <parameter name="r" value="2" dev="0.45"/>
    <parameter name="a" value="0.2" dev="0.02"/>
    <parameter name="RotBeta" value="1"/>
    <parameter name="BetaDev" value="0.2"/>
    <parameter name="InitBeta" value="0"/>
    <parameter name="tropismT" value="4"/>
    <parameter name="tropismN" value="18"/>
    <parameter name="tropismS" value="0.01"/>
    <parameter name="dx" value="0.25"/>
    <parameter name="theta" value="0" dev="0"/>
    <parameter name="rlt" value="1000000000" dev="0"/>
    <parameter name="gf" value="1"/>
    {stem}
</organ>

<organ type="stem" name="1stbranch" subType="2">
    <parameter name="lb" value="8" dev="0"/>
    <!--Basal zone [cm]-->
    <parameter name="la" value="1" dev="0"/>
    <!--Apical zone [cm];-->
    <parameter name="ln" value="10" dev="0" functiontype="1"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="lmax" value="23" dev="0"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="r" value="1"/>
    <parameter name="a" value="0.3" />
    <parameter name="RotBeta" value="1"/>
    <parameter name="BetaDev" value="0"/>
    <parameter name="InitBeta" value="0"/>
    <parameter name="tropismT" value="1"/>
    <parameter name="tropismN" value="2"/>
    <parameter name="tropismS" value="0.1"/>
    <parameter name="dx" value="0.25"/>
    <parameter name="theta" value="0" dev="0"/>
    <parameter name="rlt" value="1000000000" dev="0"/>
    <parameter name="gf" value="1"/>
</organ>
<organ type="leaf" name="lateral1" subType="1">
    <parameter name="lb" value="2" dev="0"/>
    <!--Basal zone [cm]-->
    <parameter name="la" value="3.5" dev="0"/>
    <!--Apical zone [cm];-->
    <parameter name="ln" value="0" dev="0" functiontype="3"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="lmax" value="8.5" dev="0"/>
    <!--Inter-lateral distance [cm];-->
    <parameter name="r" value="0.5" />
    <parameter name="a" value="0.016"/>
    <parameter name="tropismT" value="6"/>
    <parameter name="tropismN" value="18"/>
    <parameter name="tropismS" value="0.2"/>
    <parameter name="tropismAge" value="10"/>
    <parameter name="dx" value="1"/>
    <parameter name="dxMin" value="0.5"/>
    <parameter name="theta" value="0" dev="0"/>
    <parameter name="rlt" value="1000000000" dev="0"/>
    <parameter name="gf" value="3"/>
    <parameter name="areaMax" value="10"/>
    <parameter name="geometryN" value="100"/>
    <parameter name="parametrisationType" value="1"/>
	<parameter name="leafGeometry" phi="-3, -2.1, 0., 2.45, 3.5" x="0., 1.54, 1.7, 1.26, 0."/>
</organ>



</Plant>"""
