<?xml version="1.0"?>
<VOTABLE version="1.1" xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <INFO name="QUERY_STATUS" value="OK"/>
  <RESOURCE type="results">
    <TABLE utype="photdm:PhotometryFilter.transmissionCurve.spectrum">
    <PARAM name="FilterProfileService" value="ivo://svo/fps" ucd="meta.ref.ivorn" utype="PhotometryFilter.fpsIdentifier" datatype="char" arraysize="*"/>
    <PARAM name="filterID" value="2MASS/2MASS.Ks" ucd="meta.ref.ivoid" utype="photdm:PhotometryFilter.identifier" datatype="char" arraysize="*"/>
    <PARAM name="WavelengthUnit" value="Angstrom" ucd="meta.unit" utype="PhotometryFilter.SpectralAxis.unit" datatype="char" arraysize="*"/>
    <PARAM name="WavelengthUCD" value="em.wl" ucd="meta.ucd" utype="PhotometryFilter.SpectralAxis.UCD" datatype="char" arraysize="*"/>
    <PARAM name="Description" value="2MASS Ks" ucd="meta.note" utype="photdm:PhotometryFilter.description" datatype="char" arraysize="*"/>
    <PARAM name="PhotSystem" value="2MASS" utype="photdm:PhotometricSystem.description" datatype="char" arraysize="*">
       <DESCRIPTION>Photometric system</DESCRIPTION>
    </PARAM>
    <PARAM name="DetectorType" value="0" ucd="meta.code" utype="photdm:PhotometricSystem.detectorType" datatype="char" arraysize="*">
       <DESCRIPTION>Detector type. 0:Energy counter, 1:Photon counter.</DESCRIPTION>
    </PARAM>
    <PARAM name="Band" value="Ks" ucd="instr.bandpass" utype="photdm:PhotometryFilter.bandName" datatype="char" arraysize="*"/>
    <PARAM name="Facility" value="2MASS" ucd="instr.obsty" datatype="char" arraysize="*">
       <DESCRIPTION>Observational facility</DESCRIPTION>
    </PARAM>
    <PARAM name="ProfileReference" value="http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html#rsr" datatype="char" arraysize="*"/>
    <PARAM name="CalibrationReference" value="http://adsabs.harvard.edu/abs/2003AJ....126.1090C" datatype="char" arraysize="*"/>
    <PARAM name="Description" value="2MASS Ks" ucd="meta.note" utype="photdm:PhotometryFilter.description" datatype="char" arraysize="*"/>
    <PARAM name="components" value="Filter + Instrument + Atmosphere" datatype="char" arraysize="*">
       <DESCRIPTION>Transmission components</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthRef" value="21590" unit="Angstrom" ucd="em.wl;meta.main" utype="photdm:PhotometryFilter.spectralLocation.value" datatype="double" >
       <DESCRIPTION>Manually specified. See reference</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthMean" value="21590" unit="Angstrom" ucd="em.wl" datatype="double" >
       <DESCRIPTION>Manually specified. See reference</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthEff" value="21590" unit="Angstrom" ucd="em.wl.effective" datatype="double" >
       <DESCRIPTION>Manually specified. See reference</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthMin" value="19543.692304348" unit="Angstrom" ucd="em.wl;stat.min" utype="photdm:PhotometryFilter.bandwidth.start.value" datatype="double" >
       <DESCRIPTION>Minimum filter wavelength. Defined as the first lambda value with a transmission at least 1% of maximum transmission</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthMax" value="23552.400005256" unit="Angstrom" ucd="em.wl;stat.max" utype="photdm:PhotometryFilter.bandwidth.stop.value" datatype="double" >
       <DESCRIPTION>Maximum filter wavelength. Defined as the last lambda value with a transmission at least 1% of maximum transmission</DESCRIPTION>
    </PARAM>
    <PARAM name="WidthEff" value="2618.8695332218" unit="Angstrom" ucd="instr.bandwidth" utype="photdm:PhotometryFilter.bandwidth.extent.value" datatype="double" >
       <DESCRIPTION>Effective width. Defined as integ[x*filter(x) dx].\nEquivalent to the horizontal size of a rectangle with height equal to maximum transmission and with the same area that the one covered by the filter transmission curve.</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthCen" value="21634.039848591" unit="Angstrom" ucd="em.wl" datatype="double" >
       <DESCRIPTION>Central wavelength. Defined as the central wavelength between the two points defining FWMH</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthPivot" value="21603.094799742" unit="Angstrom" ucd="em.wl" datatype="double" >
       <DESCRIPTION>Peak wavelength. Defined as sqrt{integ[x*filter(x) dx]/integ[filter(x) dx/x]}</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthPeak" value="22630" unit="Angstrom" ucd="em.wl" datatype="double" >
       <DESCRIPTION>Peak wavelength. Defined as the lambda value with larger transmission</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthPhot" value="21557.652498909" unit="Angstrom" ucd="em.wl" datatype="double" >
       <DESCRIPTION>Photon distribution based effective wavelength. Defined as integ[x^2*filter(x)*vega(x) dx]/integ[x*filter(x)*vega(x) dx]</DESCRIPTION>
    </PARAM>
    <PARAM name="FWHM" value="2784.5527721703" unit="Angstrom" ucd="instr.bandwidth" datatype="double" >
       <DESCRIPTION>Full width at half maximum. Defined as the difference between the two wavelengths for which filter transmission is half maximum</DESCRIPTION>
    </PARAM>
    <PARAM name="Fsun" value="8.794065776272" unit="erg/cm2/s/A" ucd="phot.flux.density" datatype="double" >
       <DESCRIPTION>Sun flux</DESCRIPTION>
    </PARAM>
    <PARAM name="PhotCalID" value="2MASS/2MASS.Ks/Vega" ucd="meta.id" utype="photdm:PhotCal.identifier" datatype="char" arraysize="*"/>
    <PARAM name="MagSys" value="Vega" ucd="meta.code" utype="photdm:PhotCal.MagnitudeSystem.type" datatype="char" arraysize="*"/>
    <PARAM name="ZeroPoint" value="666.8" unit="Jy" ucd="phot.flux.density" utype="photdm:PhotCal.zeroPoint.flux.value" datatype="double" />
    <PARAM name="ZeroPointUnit" value="Jy" ucd="meta.unit" utype="photdm:PhotCal.ZeroPoint.flux.unitexpression" datatype="char" arraysize="*"/>
    <PARAM name="ZeroPointType" value="Pogson" ucd="meta.code" utype="photdm:PhotCal.ZeroPoint.type" datatype="char" arraysize="*"/>
      <FIELD name="Wavelength" utype="spec:Data.SpectralAxis.Value" ucd="em.wl" unit="Angstrom" datatype="double"/>
      <FIELD name="Transmission" utype="spec:Data.FluxAxis.Value" ucd="phys.transmission" unit="" datatype="double"/>
      <DATA>
        <TABLEDATA>
          <TR>
            <TD>19000.0</TD>
            <TD>0.0000000000</TD>
          </TR>
          <TR>
            <TD>19150.0</TD>
            <TD>0.0000000000</TD>
          </TR>
          <TR>
            <TD>19270.0</TD>
            <TD>0.0000000000</TD>
          </TR>
          <TR>
            <TD>19340.0</TD>
            <TD>0.0002000000</TD>
          </TR>
          <TR>
            <TD>19390.0</TD>
            <TD>0.0005000000</TD>
          </TR>
          <TR>
            <TD>19480.0</TD>
            <TD>0.0054000001</TD>
          </TR>
          <TR>
            <TD>19570.0</TD>
            <TD>0.0119000003</TD>
          </TR>
          <TR>
            <TD>19620.0</TD>
            <TD>0.0197000001</TD>
          </TR>
          <TR>
            <TD>19690.0</TD>
            <TD>0.0421999991</TD>
          </TR>
          <TR>
            <TD>19760.0</TD>
            <TD>0.0873000026</TD>
          </TR>
          <TR>
            <TD>19810.0</TD>
            <TD>0.1527999938</TD>
          </TR>
          <TR>
            <TD>19890.0</TD>
            <TD>0.2481999993</TD>
          </TR>
          <TR>
            <TD>19900.0</TD>
            <TD>0.1902000010</TD>
          </TR>
          <TR>
            <TD>19980.0</TD>
            <TD>0.2338999957</TD>
          </TR>
          <TR>
            <TD>20080.0</TD>
            <TD>0.2946000099</TD>
          </TR>
          <TR>
            <TD>20140.0</TD>
            <TD>0.3982000053</TD>
          </TR>
          <TR>
            <TD>20190.0</TD>
            <TD>0.3366000056</TD>
          </TR>
          <TR>
            <TD>20280.0</TD>
            <TD>0.6207000017</TD>
          </TR>
          <TR>
            <TD>20370.0</TD>
            <TD>0.7649999857</TD>
          </TR>
          <TR>
            <TD>20450.0</TD>
            <TD>0.7463999987</TD>
          </TR>
          <TR>
            <TD>20610.0</TD>
            <TD>0.6251000166</TD>
          </TR>
          <TR>
            <TD>20720.0</TD>
            <TD>0.7254999876</TD>
          </TR>
          <TR>
            <TD>20750.0</TD>
            <TD>0.6894999743</TD>
          </TR>
          <TR>
            <TD>20820.0</TD>
            <TD>0.7878999710</TD>
          </TR>
          <TR>
            <TD>20890.0</TD>
            <TD>0.8180999756</TD>
          </TR>
          <TR>
            <TD>20990.0</TD>
            <TD>0.8227999806</TD>
          </TR>
          <TR>
            <TD>21060.0</TD>
            <TD>0.8633000255</TD>
          </TR>
          <TR>
            <TD>21130.0</TD>
            <TD>0.8777999878</TD>
          </TR>
          <TR>
            <TD>21200.0</TD>
            <TD>0.8549000025</TD>
          </TR>
          <TR>
            <TD>21240.0</TD>
            <TD>0.8952999711</TD>
          </TR>
          <TR>
            <TD>21380.0</TD>
            <TD>0.9189000130</TD>
          </TR>
          <TR>
            <TD>21450.0</TD>
            <TD>0.9268000126</TD>
          </TR>
          <TR>
            <TD>21550.0</TD>
            <TD>0.9266999960</TD>
          </TR>
          <TR>
            <TD>21690.0</TD>
            <TD>0.9009000063</TD>
          </TR>
          <TR>
            <TD>21760.0</TD>
            <TD>0.9228000045</TD>
          </TR>
          <TR>
            <TD>21850.0</TD>
            <TD>0.8428000212</TD>
          </TR>
          <TR>
            <TD>21970.0</TD>
            <TD>0.9459000230</TD>
          </TR>
          <TR>
            <TD>22080.0</TD>
            <TD>0.9804000258</TD>
          </TR>
          <TR>
            <TD>22130.0</TD>
            <TD>0.9879000187</TD>
          </TR>
          <TR>
            <TD>22180.0</TD>
            <TD>0.9847999811</TD>
          </TR>
          <TR>
            <TD>22320.0</TD>
            <TD>0.9646999836</TD>
          </TR>
          <TR>
            <TD>22370.0</TD>
            <TD>0.9815999866</TD>
          </TR>
          <TR>
            <TD>22480.0</TD>
            <TD>0.9833999872</TD>
          </TR>
          <TR>
            <TD>22560.0</TD>
            <TD>0.9613000154</TD>
          </TR>
          <TR>
            <TD>22600.0</TD>
            <TD>0.9792000055</TD>
          </TR>
          <TR>
            <TD>22630.0</TD>
            <TD>1.0000000000</TD>
          </TR>
          <TR>
            <TD>22650.0</TD>
            <TD>0.9631999731</TD>
          </TR>
          <TR>
            <TD>22700.0</TD>
            <TD>0.9811999798</TD>
          </TR>
          <TR>
            <TD>22720.0</TD>
            <TD>0.9681000113</TD>
          </TR>
          <TR>
            <TD>22760.0</TD>
            <TD>0.9108999968</TD>
          </TR>
          <TR>
            <TD>22770.0</TD>
            <TD>0.9821000099</TD>
          </TR>
          <TR>
            <TD>22810.0</TD>
            <TD>0.8895999789</TD>
          </TR>
          <TR>
            <TD>22840.0</TD>
            <TD>0.8917999864</TD>
          </TR>
          <TR>
            <TD>22860.0</TD>
            <TD>0.9423999786</TD>
          </TR>
          <TR>
            <TD>22910.0</TD>
            <TD>0.8403999805</TD>
          </TR>
          <TR>
            <TD>22930.0</TD>
            <TD>0.8041999936</TD>
          </TR>
          <TR>
            <TD>22950.0</TD>
            <TD>0.7077000141</TD>
          </TR>
          <TR>
            <TD>22970.0</TD>
            <TD>0.6575999856</TD>
          </TR>
          <TR>
            <TD>22990.0</TD>
            <TD>0.5606999993</TD>
          </TR>
          <TR>
            <TD>23060.0</TD>
            <TD>0.4436999857</TD>
          </TR>
          <TR>
            <TD>23110.0</TD>
            <TD>0.3481999934</TD>
          </TR>
          <TR>
            <TD>23160.0</TD>
            <TD>0.2301999927</TD>
          </TR>
          <TR>
            <TD>23200.0</TD>
            <TD>0.1625999957</TD>
          </TR>
          <TR>
            <TD>23250.0</TD>
            <TD>0.1360000074</TD>
          </TR>
          <TR>
            <TD>23280.0</TD>
            <TD>0.0921000019</TD>
          </TR>
          <TR>
            <TD>23350.0</TD>
            <TD>0.0623999983</TD>
          </TR>
          <TR>
            <TD>23390.0</TD>
            <TD>0.0430999994</TD>
          </TR>
          <TR>
            <TD>23440.0</TD>
            <TD>0.0340000018</TD>
          </TR>
          <TR>
            <TD>23460.0</TD>
            <TD>0.0309999995</TD>
          </TR>
          <TR>
            <TD>23520.0</TD>
            <TD>0.0118000004</TD>
          </TR>
          <TR>
            <TD>23610.0</TD>
            <TD>0.0068000001</TD>
          </TR>
          <TR>
            <TD>23630.0</TD>
            <TD>0.0007000000</TD>
          </TR>
          <TR>
            <TD>23700.0</TD>
            <TD>0.0030000000</TD>
          </TR>
          <TR>
            <TD>23750.0</TD>
            <TD>0.0020999999</TD>
          </TR>
          <TR>
            <TD>23840.0</TD>
            <TD>0.0004000000</TD>
          </TR>
          <TR>
            <TD>23990.0</TD>
            <TD>0.0000000000</TD>
          </TR>
        </TABLEDATA>
      </DATA>
    </TABLE>
  </RESOURCE>
</VOTABLE>
