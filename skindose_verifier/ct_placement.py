"""Estimated placement of a Slicer surface using a CT-derived target point."""
import numpy as np
from pydicom import dcmread
from trimesh import load_mesh

from .geometry import beam_from_event, table_shift


def _coded_values(root, concept):
    values = set()
    for item in getattr(root, "ContentSequence", ()):
        name = getattr(item, "ConceptNameCodeSequence", ())
        value = getattr(item, "ConceptCodeSequence", ())
        if name and value and str(name[0].CodeMeaning).strip().casefold() == concept:
            values.add(str(value[0].CodeMeaning).strip().casefold())
        values.update(_coded_values(item, concept))
    return values


def estimate_placement(path, digest, metadata, dicom, events, equipment):
    if metadata.get("mesh_file") != path.name or metadata.get("mesh_sha256") != digest:
        raise ValueError("El JSON de Slicer no corresponde al STL seleccionado.")
    if metadata.get("coordinate_system") != "RAS" or metadata.get("length_unit") != "mm":
        raise ValueError("La colocación automática requiere un STL documentado en RAS y mm.")
    sr = dcmread(dicom)
    direction = _coded_values(sr, "patient table relationship")
    posture = _coded_values(sr, "patient orientation modifier")
    if direction not in ({"headfirst"}, {"feetfirst"}) or posture not in ({"supine"}, {"prone"}):
        raise ValueError("El RDSR no identifica una postura supina o prona uniforme y el sentido de entrada.")
    regions = _coded_values(sr, "target region")
    anchors = metadata.get("anatomical_anchors_ras_mm", {})
    if len(regions) != 1 or not isinstance(anchors, dict) or not regions.issubset(anchors):
        raise ValueError("Falta un punto anatómico derivado del TC para la región diana del RDSR.")
    region = next(iter(regions))
    anchor = np.asarray(anchors[region], dtype=float)
    mesh = load_mesh(path, process=False)
    bounds = mesh.bounds
    if anchor.shape != (3,) or not np.isfinite(anchor).all() or bounds is None:
        raise ValueError("Punto anatómico o superficie inválidos.")
    if np.any(anchor < bounds[0]) or np.any(anchor > bounds[1]):
        raise ValueError("El punto de la región diana queda fuera de los límites del STL.")
    # RAS to left, cranial, anterior for head-first supine.
    matrix = np.array([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])
    if direction == {"feetfirst"}:
        matrix = np.diag([-1, -1, 1]) @ matrix
    if posture == {"prone"}:
        matrix = np.diag([-1, 1, -1]) @ matrix
    if not events:
        raise ValueError("No hay eventos para establecer la posición de referencia.")
    event = events[0]
    shift = table_shift(event, equipment)
    beam = beam_from_event(event, equipment)
    geometry = equipment.geometry_for(event.plane)
    support_z = float(geometry["table_plane_z_mm"])
    local_anchor = matrix @ anchor
    min_z = float(np.min(mesh.vertices @ matrix[2]))
    offset_z = support_z - min_z
    room_anchor_z = local_anchor[2] + offset_z + shift[2]
    iso = beam.source_mm + float(event.sod_mm) * beam.direction
    if abs(beam.direction[2]) < 0.1:
        raise ValueError("El haz de referencia es casi paralelo a la mesa; el apoyo no fija el anclaje.")
    depth = (room_anchor_z - beam.source_mm[2]) / beam.direction[2]
    if not 0 < depth < float(event.sid_mm):
        raise ValueError("El anclaje estimado queda fuera del tramo entre foco y detector.")
    target = beam.source_mm + depth * beam.direction
    offset = target - shift - local_anchor
    return dict(mesh_unit="mm", mesh_axes_matrix=matrix.tolist(),
                mesh_origin_mm=offset.tolist(), origin_on_table_mm=[0, 0, 0],
                tilt_x_deg=0, tilt_y_deg=0), dict(
        target_region=region, anchor_ras_mm=anchor.tolist(),
        reference_event_index=1, reference_table_shift_mm=shift.tolist(),
        reference_isocenter_mm=iso.tolist(), posture=next(iter(posture)),
        table_relationship=next(iter(direction)),
        assumptions=["superficie apoyada en el plano de mesa configurado",
                     "punto diana en el eje central del primer evento",
                     "paciente fijo respecto a la mesa"])
