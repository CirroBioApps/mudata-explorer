from mudata_explorer.base import transform as transform_dir
from mudata_explorer.base.transform import Transform
from typing import Dict


def all_transforms() -> Dict[str, Transform]:
    return {
        getattr(transform_dir, transform).id: getattr(transform_dir, transform)
        for transform in dir(transform_dir)
        if not transform.startswith("__")
        and hasattr(getattr(transform_dir, transform), "id")
    }


def get_transform(id: str) -> Transform:
    transforms = all_transforms()
    if id not in transforms:
        raise ValueError(f"Transform with id '{id}' not found.")
    return transforms[id]
