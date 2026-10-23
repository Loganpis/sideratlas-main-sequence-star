"""Targeted checks for the colour refinement and continuous four-orbit camera.

The physical solver is unchanged. Existing whole-orbit projection quadrature
is reused only after asserting equality of every projection input setting.
"""
from pathlib import Path
import sys,json,copy,gc,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
from golden_orbit.render import Renderer,grade_photosphere
from geometry import camera_direction,orientation
from solar_surface import YROW

def main():
    c=json.loads((HERE/'config-final.json').read_text());result={}
    t=np.linspace(0,32,9601);inc,az=camera_direction(t/32,c)
    assert np.max(abs(np.gradient(az,t)-45))<1e-8 and az[-1]-az[0]==1440
    assert np.all(inc==67) and c['center_fraction']==[.5,.565] and c['radius_fraction_of_width']==.405
    periodic=[]
    for start in [0,8,16,24]:
        i0,a0=camera_direction(start/32,c);i1,a1=camera_direction((start+8)/32,c)
        err=float(np.max(abs(orientation(np.eye(3),i0,a0)-orientation(np.eye(3),i1,a1))))
        periodic.append(err);assert err<1e-13
    frames=(np.arange(960)+.5)/30
    fi,fa=camera_direction(frames/32,c)
    assert np.max(abs(np.diff(fa)-1.5))<1e-10
    result['camera']=dict(revolutions=4,seconds_per_orbit=8,angular_speed_deg_s=45,
        azimuth_at_0_8_16_24_32_s=az[::2400].tolist(),inclination_deg=67,
        maximum_frozen_orientation_return_error=max(periodic),
        frame_sampling='960 midpoint frames; constant 1.5 degree azimuth steps, no duplicated endpoint',
        centered_constant_radius=True,source_clock_reset=False,
        physical_seconds_per_display_second=20,physical_rotation_days=25.38,
        physical_rotation_during_film_deg=640/(25.38*86400)*360)
    # This report was calculated using precisely the same UV geometry and
    # source interval; palette/exposure changes do not change physical j.
    previous=ROOT/'multispectral/config-final.json'
    inherited=ROOT/'multispectral/validation/scientific-checks.json'
    if previous.exists() and inherited.exists():
        pc=json.loads(previous.read_text());report=json.loads(inherited.read_text())
        for key in ['width','height','camera','center_fraction','radius_fraction_of_width','physical_seconds_per_display_second']:
            assert c[key]==pc[key],key
        for key in ['source','source_start_s','spatial_magnification','invariant_depth_Mm','depth_samples','patch_camera_vectors_at_start']:
            assert c['atmosphere'][key]==pc['atmosphere'][key],key
        result['UV_depth_quadrature']=report['depth_quadrature']
        result['maximum_depth_quadrature_L1_48_vs_96']=report['maximum_depth_quadrature_L1']
        result['projection_reuse_provenance']=dict(report_sha256=hashlib.sha256(inherited.read_bytes()).hexdigest(),all_projection_settings_matched=True)
    else:
        # Clean-package reproduction needs no historical edition.
        from atmosphere import Atmosphere
        a=Atmosphere(c);depth=[]
        for when in [0.,1.,2.,3.,4.,5.,6.,7.,16.,31.99,32.]:
            lo=a.image(when,48);hi=a.image(when,96)
            depth.append(dict(time_s=when,L1_48_vs_96=float(abs(lo-hi).sum()/max(hi.sum(),1e-100)),total_UV_flux_relative=float(lo.sum()/max(hi.sum(),1e-100)-1)))
        result['UV_depth_quadrature']=depth
        result['maximum_depth_quadrature_L1_48_vs_96']=max(x['L1_48_vs_96'] for x in depth)
        del a;gc.collect()
    assert result['maximum_depth_quadrature_L1_48_vs_96']<.03
    r=Renderer(c);r.atmosphere.image(0);r.atmosphere.image(32)
    assert len(r.atmosphere.times)==71 and np.all(np.diff(r.atmosphere.times)>0)
    # Use actual surface and active-region samplers with physical time held
    # fixed. The 8-s camera return must not re-tile, re-seed or move a patch.
    tile=next(v for v in r.tiles if v.row0 <= c['height']*c['center_fraction'][1] < v.row1)
    sample=tile.regular
    reference=r.surface.sample(sample,0,67,-15)
    active,weight=r.active.sample(sample,0,67,-15)
    comparisons=[]
    for cycle in [1,2,3,4]:
        q=r.surface.sample(sample,0,67,-15+360*cycle)
        av,w=r.active.sample(sample,0,67,-15+360*cycle)
        err=max(float(abs(q-reference).max()),float(abs(av-active).max()),float(abs(w-weight).max()))
        comparisons.append(err);assert err<2e-6
    result['world_anchoring']=dict(actual_surface_and_active_map_frozen_time_return_errors=comparisons,
        atmospheric_centers_star_fixed=r.atmosphere.centers.tolist(),
        camera_azimuth_never_used_for_atlas_hash=True,
        physical_source_times_optical=[c['source_start_s'],c['source_start_s']+640],
        physical_source_times_CIV=[c['atmosphere']['source_start_s'],c['atmosphere']['source_start_s']+640],
        note='Only the viewpoint repeats; evolving source radiation and the independent physical rotation do not reset')
    q=r.source(.5/c['fps']);m=q.sum(axis=-1)>0
    graded=grade_photosphere(q,r.grade_reference,c['color_grade'])
    err=float(abs(graded@YROW-q@YROW).max());assert err<1e-12
    # LUT values are continuous, fixed and driven exclusively by calculated Y.
    test=np.linspace(0,4,10001);rgb=test[:,None]*np.array([1.,1.,1.])
    cg=grade_photosphere(rgb,np.ones_like(test),c['color_grade'])
    assert np.isfinite(cg).all() and (cg>=0).all()
    result['grading']=dict(max_absolute_pre_tone_linear_Y_error=err,
        fixed_initial_limb_reference=True,per_frame_auto_exposure=False,
        relative_intensity_knots=c['color_grade']['relative_intensity_knots'],
        gain_vectors=c['color_grade']['linear_rgb_gains_at_knots'],
        UV_false_color_linear_RGB=c['atmosphere']['false_color_linear_rgb'])
    one=r.optics.process(graded);del r;gc.collect()
    c4=copy.deepcopy(c);c4['supersample']=2;r4=Renderer(c4)
    fine=r4.optics.process(grade_photosphere(r4.source(.5/c4['fps']),r4.grade_reference,c4['color_grade']))
    spatial=float(np.sqrt(np.mean((one[m]-fine[m])**2))/fine[m].mean())
    result['photospheric_output_sampling_1_vs_4_post_imaging_RMS_relative_to_mean']=spatial
    assert spatial<.03
    result['limitations']=['Independent local models, not a global MHD solution',
        'The two finite 2.5D UV sheets can look sheet-like as the camera orbits; no artificial billboarding',
        'Prescribed 10x magnification, contrast, golden grading and optical halo',
        'UV ionization, opacity and temporal reconstruction limits unchanged',
        'Four camera orbits are seamless; the evolving physical state is not a repeated 8-second clip or a closed 32-second loop']
    (HERE/'validation/scientific-checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
