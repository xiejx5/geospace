"""Bundle tested GDAL wheel payloads into GeoSpace wheels."""

from __future__ import annotations

import base64
import copy
import csv
import hashlib
import io
from pathlib import Path
import shutil
import zipfile


BUILD = Path("build")
GDAL = Path("gdal-wheels")
DIST = Path("dist")


def hash_value(data):
    digest = hashlib.sha256(data).digest()
    return "sha256=" + base64.urlsafe_b64encode(digest).rstrip(b"=").decode()


def write_wheel(path, files, record):
    rows = []

    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as wheel:
        for name, (info, data) in files.items():
            info = copy.copy(info)
            info.filename = name
            wheel.writestr(info, data)
            rows.append((name, hash_value(data), str(len(data))))

        text = io.StringIO()
        writer = csv.writer(text, lineterminator="\n")
        writer.writerows(rows)
        writer.writerow((record, "", ""))
        wheel.writestr(record, text.getvalue())


def main():
    wheels = list(BUILD.glob("geospace-*-py3-none-any.whl"))
    if len(wheels) != 1:
        raise RuntimeError(f"Expected one GeoSpace wheel, found {wheels}")

    with zipfile.ZipFile(wheels[0]) as wheel:
        geospace = {
            info.filename: (info, wheel.read(info))
            for info in wheel.infolist()
        }

    infos = {
        name.split("/", 1)[0]
        for name in geospace
        if ".dist-info/" in name
    }
    if len(infos) != 1:
        raise RuntimeError(f"Unexpected GeoSpace metadata: {infos}")

    info_dir = infos.pop()
    version = info_dir.removeprefix("geospace-").removesuffix(".dist-info")
    wheel_meta = f"{info_dir}/WHEEL"
    metadata = f"{info_dir}/METADATA"
    record = f"{info_dir}/RECORD"

    shutil.rmtree(DIST, ignore_errors=True)
    DIST.mkdir()

    for source in sorted(GDAL.glob("gdal-*.whl")):
        files = {}

        with zipfile.ZipFile(source) as wheel:
            gdal_info = next(
                name.split("/", 1)[0]
                for name in wheel.namelist()
                if ".dist-info/WHEEL" in name
            )
            gdal_data = gdal_info.removesuffix(".dist-info") + ".data/"

            tags = [
                line
                for line in wheel.read(f"{gdal_info}/WHEEL").decode().splitlines()
                if line.startswith("Tag: ")
            ]

            for item in wheel.infolist():
                name = item.filename

                if name.startswith(gdal_info + "/"):
                    continue

                if name.startswith(gdal_data):
                    name = f"geospace-{version}.data/" + name[len(gdal_data):]

                files[name] = (item, wheel.read(item))

        for name, value in geospace.items():
            if name != record:
                files[name] = value

        lines = files[wheel_meta][1].decode().splitlines()
        lines = [
            line
            for line in lines
            if not line.startswith(("Root-Is-Purelib:", "Tag: "))
        ]
        lines += ["Root-Is-Purelib: false", *tags, ""]
        files[wheel_meta] = (files[wheel_meta][0], "\n".join(lines).encode())

        required = {
            "geospace/__init__.py",
            "osgeo/__init__.py",
            "osgeo_utils/__init__.py",
            "osgeo/data/gdal/gdalvrt.xsd",
            "osgeo/data/proj/proj.db",
        }
        missing = required - files.keys()
        if missing:
            raise RuntimeError(f"{source.name}: missing {sorted(missing)}")

        if b"requires-dist: gdal" in files[metadata][1].lower():
            raise RuntimeError("GeoSpace still declares a gdal dependency")

        _, python, abi, platform = source.name[:-4].rsplit("-", 3)
        target = DIST / f"geospace-{version}-{python}-{abi}-{platform}.whl"

        write_wheel(target, files, record)
        print(target)


if __name__ == "__main__":
    main()
