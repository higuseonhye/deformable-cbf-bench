"""Isaac Lab scene sketch. Not imported by tests and not executed in this repo.

Isaac Lab 3 moved deformable materials into the backend packages. The
control code never reads this file. It consumes an (N, 3) vertex array
through deformable_cbf.tracker.VertexEllipsoidTracker.

Launching the function below requires a local Isaac Lab 3 install with a
PhysX GPU backend. That install is not part of requirements.txt.
"""

from __future__ import annotations


def build_scene_cfg():
    """Return an InteractiveSceneCfg, or raise ImportError with the reason."""
    try:
        import isaaclab.sim as sim_utils
        from isaaclab.assets import AssetBaseCfg, DeformableObjectCfg
        from isaaclab.scene import InteractiveSceneCfg
        from isaaclab.utils import configclass
        from isaaclab_physx.sim import PhysxDeformableBodyMaterialCfg
    except ImportError as exc:
        raise ImportError(
            "Isaac Lab 3 is not installed. The kinematic filter does not need it. "
            "VertexEllipsoidTracker accepts the nodal positions once a simulator is running."
        ) from exc

    @configclass
    class SurgicalDeformableSceneCfg(InteractiveSceneCfg):
        ground = AssetBaseCfg(prim_path="/World/defaultGroundPlane", spawn=sim_utils.GroundPlaneCfg())
        light = AssetBaseCfg(
            prim_path="/World/Light",
            spawn=sim_utils.DomeLightCfg(intensity=3000.0, color=(0.95, 0.95, 1.0)),
        )
        tissue = DeformableObjectCfg(
            prim_path="{ENV_REGEX_NS}/Tissue",
            spawn=sim_utils.MeshCuboidCfg(
                size=(0.08, 0.08, 0.03),
                volume_deformable_props=[],
                visual_material=sim_utils.PreviewSurfaceCfg(diffuse_color=(0.7, 0.15, 0.15)),
                physics_material=PhysxDeformableBodyMaterialCfg(
                    youngs_modulus=5.0e4,
                    poissons_ratio=0.45,
                    density=1000.0,
                ),
            ),
            init_state=DeformableObjectCfg.InitialStateCfg(pos=(0.0, 0.0, 0.015)),
        )

    return SurgicalDeformableSceneCfg
