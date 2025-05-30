<?xml version="1.0"?>
<VOTABLE version="1.1" xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <INFO name="QUERY_STATUS" value="OK"/>
  <RESOURCE type="results">
    <TABLE utype="photdm:PhotometryFilter.transmissionCurve.spectrum">
    <PARAM name="FilterProfileService" value="ivo://svo/fps" ucd="meta.ref.ivorn" utype="PhotometryFilter.fpsIdentifier" datatype="char" arraysize="*"/>
    <PARAM name="filterID" value="SLOAN/SDSS.g" ucd="meta.ref.ivoid" utype="photdm:PhotometryFilter.identifier" datatype="char" arraysize="*"/>
    <PARAM name="WavelengthUnit" value="Angstrom" ucd="meta.unit" utype="PhotometryFilter.SpectralAxis.unit" datatype="char" arraysize="*"/>
    <PARAM name="WavelengthUCD" value="em.wl" ucd="meta.ucd" utype="PhotometryFilter.SpectralAxis.UCD" datatype="char" arraysize="*"/>
    <PARAM name="Description" value="SDSS g full transmission" ucd="meta.note" utype="photdm:PhotometryFilter.description" datatype="char" arraysize="*"/>
    <PARAM name="PhotSystem" value="SDSS" utype="photdm:PhotometricSystem.description" datatype="char" arraysize="*">
       <DESCRIPTION>Photometric system</DESCRIPTION>
    </PARAM>
    <PARAM name="DetectorType" value="1" ucd="meta.code" utype="photdm:PhotometricSystem.detectorType" datatype="char" arraysize="*">
       <DESCRIPTION>Detector type. 0:Energy counter, 1:Photon counter.</DESCRIPTION>
    </PARAM>
    <PARAM name="Band" value="g" ucd="instr.bandpass" utype="photdm:PhotometryFilter.bandName" datatype="char" arraysize="*"/>
    <PARAM name="Facility" value="SLOAN" ucd="instr.obsty" datatype="char" arraysize="*">
       <DESCRIPTION>Observational facility</DESCRIPTION>
    </PARAM>
    <PARAM name="ProfileReference" value="http://www.sdss.org/dr7/instruments/imager/index.html" datatype="char" arraysize="*"/>
    <PARAM name="CalibrationReference" value="http://www.sdss.org/DR2/algorithms/fluxcal.html" datatype="char" arraysize="*"/>
    <PARAM name="Description" value="SDSS g full transmission" ucd="meta.note" utype="photdm:PhotometryFilter.description" datatype="char" arraysize="*"/>
    <PARAM name="components" value="Filter + Instrument + Atmosphere" datatype="char" arraysize="*">
       <DESCRIPTION>Transmission components</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthRef" value="4702.4953002767" unit="Angstrom" ucd="em.wl;meta.main" utype="photdm:PhotometryFilter.spectralLocation.value" datatype="double" >
       <DESCRIPTION>Reference wavelength. Defined as the same than the pivot wavelength.</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthMean" value="4750.8231056844" unit="Angstrom" ucd="em.wl" datatype="double" >
       <DESCRIPTION>Mean wavelength. Defined as integ[x*filter(x) dx]/integ[filter(x) dx]</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthEff" value="4671.7822137652" unit="Angstrom" ucd="em.wl.effective" datatype="double" >
       <DESCRIPTION>Effective wavelength. Defined as integ[x*filter(x)*vega(x) dx]/integ[filter(x)*vega(x) dx]</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthMin" value="3797.6384743979" unit="Angstrom" ucd="em.wl;stat.min" utype="photdm:PhotometryFilter.bandwidth.start.value" datatype="double" >
       <DESCRIPTION>Minimum filter wavelength. Defined as the first lambda value with a transmission at least 1% of maximum transmission</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthMax" value="5553.0413712781" unit="Angstrom" ucd="em.wl;stat.max" utype="photdm:PhotometryFilter.bandwidth.stop.value" datatype="double" >
       <DESCRIPTION>Maximum filter wavelength. Defined as the last lambda value with a transmission at least 1% of maximum transmission</DESCRIPTION>
    </PARAM>
    <PARAM name="WidthEff" value="1064.6831251068" unit="Angstrom" ucd="instr.bandwidth" utype="photdm:PhotometryFilter.bandwidth.extent.value" datatype="double" >
       <DESCRIPTION>Effective width. Defined as integ[x*filter(x) dx].\nEquivalent to the horizontal size of a rectangle with height equal to maximum transmission and with the same area that the one covered by the filter transmission curve.</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthCen" value="4747.3276268194" unit="Angstrom" ucd="em.wl" datatype="double" >
       <DESCRIPTION>Central wavelength. Defined as the central wavelength between the two points defining FWMH</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthPivot" value="4702.4953002767" unit="Angstrom" ucd="em.wl" datatype="double" >
       <DESCRIPTION>Peak wavelength. Defined as sqrt{integ[x*filter(x) dx]/integ[filter(x) dx/x]}</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthPeak" value="5180" unit="Angstrom" ucd="em.wl" datatype="double" >
       <DESCRIPTION>Peak wavelength. Defined as the lambda value with larger transmission</DESCRIPTION>
    </PARAM>
    <PARAM name="WavelengthPhot" value="4703.954046074" unit="Angstrom" ucd="em.wl" datatype="double" >
       <DESCRIPTION>Photon distribution based effective wavelength. Defined as integ[x^2*filter(x)*vega(x) dx]/integ[x*filter(x)*vega(x) dx]</DESCRIPTION>
    </PARAM>
    <PARAM name="FWHM" value="1175.631602652" unit="Angstrom" ucd="instr.bandwidth" datatype="double" >
       <DESCRIPTION>Full width at half maximum. Defined as the difference between the two wavelengths for which filter transmission is half maximum</DESCRIPTION>
    </PARAM>
    <PARAM name="Fsun" value="187.21772917722" unit="erg/cm2/s/A" ucd="phot.flux.density" datatype="double" >
       <DESCRIPTION>Sun flux</DESCRIPTION>
    </PARAM>
    <PARAM name="PhotCalID" value="SLOAN/SDSS.g/Vega" ucd="meta.id" utype="photdm:PhotCal.identifier" datatype="char" arraysize="*"/>
    <PARAM name="MagSys" value="Vega" ucd="meta.code" utype="photdm:PhotCal.MagnitudeSystem.type" datatype="char" arraysize="*"/>
    <PARAM name="ZeroPoint" value="4023.5732569791" unit="Jy" ucd="phot.flux.density" utype="photdm:PhotCal.zeroPoint.flux.value" datatype="double" />
    <PARAM name="ZeroPointUnit" value="Jy" ucd="meta.unit" utype="photdm:PhotCal.ZeroPoint.flux.unitexpression" datatype="char" arraysize="*"/>
    <PARAM name="ZeroPointType" value="Pogson" ucd="meta.code" utype="photdm:PhotCal.ZeroPoint.type" datatype="char" arraysize="*"/>
      <FIELD name="Wavelength" utype="spec:Data.SpectralAxis.Value" ucd="em.wl" unit="Angstrom" datatype="double"/>
      <FIELD name="Transmission" utype="spec:Data.FluxAxis.Value" ucd="phys.transmission" unit="" datatype="double"/>
      <DATA>
        <TABLEDATA>
          <TR>
            <TD>3630.0</TD>
            <TD>0.0000000000</TD>
          </TR>
          <TR>
            <TD>3655.0</TD>
            <TD>0.0003000000</TD>
          </TR>
          <TR>
            <TD>3680.0</TD>
            <TD>0.0008000000</TD>
          </TR>
          <TR>
            <TD>3705.0</TD>
            <TD>0.0013000000</TD>
          </TR>
          <TR>
            <TD>3730.0</TD>
            <TD>0.0019000000</TD>
          </TR>
          <TR>
            <TD>3755.0</TD>
            <TD>0.0024000000</TD>
          </TR>
          <TR>
            <TD>3780.0</TD>
            <TD>0.0034000000</TD>
          </TR>
          <TR>
            <TD>3805.0</TD>
            <TD>0.0055000000</TD>
          </TR>
          <TR>
            <TD>3830.0</TD>
            <TD>0.0103000000</TD>
          </TR>
          <TR>
            <TD>3855.0</TD>
            <TD>0.0194000000</TD>
          </TR>
          <TR>
            <TD>3880.0</TD>
            <TD>0.0326000000</TD>
          </TR>
          <TR>
            <TD>3905.0</TD>
            <TD>0.0492000000</TD>
          </TR>
          <TR>
            <TD>3930.0</TD>
            <TD>0.0686000000</TD>
          </TR>
          <TR>
            <TD>3955.0</TD>
            <TD>0.0900000000</TD>
          </TR>
          <TR>
            <TD>3980.0</TD>
            <TD>0.1123000000</TD>
          </TR>
          <TR>
            <TD>4005.0</TD>
            <TD>0.1342000000</TD>
          </TR>
          <TR>
            <TD>4030.0</TD>
            <TD>0.1545000000</TD>
          </TR>
          <TR>
            <TD>4055.0</TD>
            <TD>0.1722000000</TD>
          </TR>
          <TR>
            <TD>4080.0</TD>
            <TD>0.1873000000</TD>
          </TR>
          <TR>
            <TD>4105.0</TD>
            <TD>0.2003000000</TD>
          </TR>
          <TR>
            <TD>4130.0</TD>
            <TD>0.2116000000</TD>
          </TR>
          <TR>
            <TD>4155.0</TD>
            <TD>0.2214000000</TD>
          </TR>
          <TR>
            <TD>4180.0</TD>
            <TD>0.2301000000</TD>
          </TR>
          <TR>
            <TD>4205.0</TD>
            <TD>0.2378000000</TD>
          </TR>
          <TR>
            <TD>4230.0</TD>
            <TD>0.2448000000</TD>
          </TR>
          <TR>
            <TD>4255.0</TD>
            <TD>0.2513000000</TD>
          </TR>
          <TR>
            <TD>4280.0</TD>
            <TD>0.2574000000</TD>
          </TR>
          <TR>
            <TD>4305.0</TD>
            <TD>0.2633000000</TD>
          </TR>
          <TR>
            <TD>4330.0</TD>
            <TD>0.2691000000</TD>
          </TR>
          <TR>
            <TD>4355.0</TD>
            <TD>0.2747000000</TD>
          </TR>
          <TR>
            <TD>4380.0</TD>
            <TD>0.2801000000</TD>
          </TR>
          <TR>
            <TD>4405.0</TD>
            <TD>0.2852000000</TD>
          </TR>
          <TR>
            <TD>4430.0</TD>
            <TD>0.2899000000</TD>
          </TR>
          <TR>
            <TD>4455.0</TD>
            <TD>0.2940000000</TD>
          </TR>
          <TR>
            <TD>4480.0</TD>
            <TD>0.2979000000</TD>
          </TR>
          <TR>
            <TD>4505.0</TD>
            <TD>0.3016000000</TD>
          </TR>
          <TR>
            <TD>4530.0</TD>
            <TD>0.3055000000</TD>
          </TR>
          <TR>
            <TD>4555.0</TD>
            <TD>0.3097000000</TD>
          </TR>
          <TR>
            <TD>4580.0</TD>
            <TD>0.3141000000</TD>
          </TR>
          <TR>
            <TD>4605.0</TD>
            <TD>0.3184000000</TD>
          </TR>
          <TR>
            <TD>4630.0</TD>
            <TD>0.3224000000</TD>
          </TR>
          <TR>
            <TD>4655.0</TD>
            <TD>0.3257000000</TD>
          </TR>
          <TR>
            <TD>4680.0</TD>
            <TD>0.3284000000</TD>
          </TR>
          <TR>
            <TD>4705.0</TD>
            <TD>0.3307000000</TD>
          </TR>
          <TR>
            <TD>4730.0</TD>
            <TD>0.3327000000</TD>
          </TR>
          <TR>
            <TD>4755.0</TD>
            <TD>0.3346000000</TD>
          </TR>
          <TR>
            <TD>4780.0</TD>
            <TD>0.3364000000</TD>
          </TR>
          <TR>
            <TD>4805.0</TD>
            <TD>0.3383000000</TD>
          </TR>
          <TR>
            <TD>4830.0</TD>
            <TD>0.3403000000</TD>
          </TR>
          <TR>
            <TD>4855.0</TD>
            <TD>0.3425000000</TD>
          </TR>
          <TR>
            <TD>4880.0</TD>
            <TD>0.3448000000</TD>
          </TR>
          <TR>
            <TD>4905.0</TD>
            <TD>0.3472000000</TD>
          </TR>
          <TR>
            <TD>4930.0</TD>
            <TD>0.3495000000</TD>
          </TR>
          <TR>
            <TD>4955.0</TD>
            <TD>0.3519000000</TD>
          </TR>
          <TR>
            <TD>4980.0</TD>
            <TD>0.3541000000</TD>
          </TR>
          <TR>
            <TD>5005.0</TD>
            <TD>0.3562000000</TD>
          </TR>
          <TR>
            <TD>5030.0</TD>
            <TD>0.3581000000</TD>
          </TR>
          <TR>
            <TD>5055.0</TD>
            <TD>0.3597000000</TD>
          </TR>
          <TR>
            <TD>5080.0</TD>
            <TD>0.3609000000</TD>
          </TR>
          <TR>
            <TD>5105.0</TD>
            <TD>0.3613000000</TD>
          </TR>
          <TR>
            <TD>5130.0</TD>
            <TD>0.3609000000</TD>
          </TR>
          <TR>
            <TD>5155.0</TD>
            <TD>0.3595000000</TD>
          </TR>
          <TR>
            <TD>5180.0</TD>
            <TD>0.3581000000</TD>
          </TR>
          <TR>
            <TD>5205.0</TD>
            <TD>0.3558000000</TD>
          </TR>
          <TR>
            <TD>5230.0</TD>
            <TD>0.3452000000</TD>
          </TR>
          <TR>
            <TD>5255.0</TD>
            <TD>0.3194000000</TD>
          </TR>
          <TR>
            <TD>5280.0</TD>
            <TD>0.2807000000</TD>
          </TR>
          <TR>
            <TD>5305.0</TD>
            <TD>0.2339000000</TD>
          </TR>
          <TR>
            <TD>5330.0</TD>
            <TD>0.1839000000</TD>
          </TR>
          <TR>
            <TD>5355.0</TD>
            <TD>0.1352000000</TD>
          </TR>
          <TR>
            <TD>5380.0</TD>
            <TD>0.0911000000</TD>
          </TR>
          <TR>
            <TD>5405.0</TD>
            <TD>0.0548000000</TD>
          </TR>
          <TR>
            <TD>5430.0</TD>
            <TD>0.0295000000</TD>
          </TR>
          <TR>
            <TD>5455.0</TD>
            <TD>0.0166000000</TD>
          </TR>
          <TR>
            <TD>5480.0</TD>
            <TD>0.0112000000</TD>
          </TR>
          <TR>
            <TD>5505.0</TD>
            <TD>0.0077000000</TD>
          </TR>
          <TR>
            <TD>5530.0</TD>
            <TD>0.0050000000</TD>
          </TR>
          <TR>
            <TD>5555.0</TD>
            <TD>0.0032000000</TD>
          </TR>
          <TR>
            <TD>5580.0</TD>
            <TD>0.0021000000</TD>
          </TR>
          <TR>
            <TD>5605.0</TD>
            <TD>0.0015000000</TD>
          </TR>
          <TR>
            <TD>5630.0</TD>
            <TD>0.0012000000</TD>
          </TR>
          <TR>
            <TD>5655.0</TD>
            <TD>0.0010000000</TD>
          </TR>
          <TR>
            <TD>5680.0</TD>
            <TD>0.0009000000</TD>
          </TR>
          <TR>
            <TD>5705.0</TD>
            <TD>0.0008000000</TD>
          </TR>
          <TR>
            <TD>5730.0</TD>
            <TD>0.0006000000</TD>
          </TR>
          <TR>
            <TD>5755.0</TD>
            <TD>0.0005000000</TD>
          </TR>
          <TR>
            <TD>5780.0</TD>
            <TD>0.0003000000</TD>
          </TR>
          <TR>
            <TD>5805.0</TD>
            <TD>0.0001000000</TD>
          </TR>
          <TR>
            <TD>5830.0</TD>
            <TD>0.0000000000</TD>
          </TR>
        </TABLEDATA>
      </DATA>
    </TABLE>
  </RESOURCE>
</VOTABLE>
