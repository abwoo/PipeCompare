"""Metadata-only OSF discovery stays inside public API URLs and pinned hashes."""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks"))
from discover_holdout_osf_metadata import discover,endpoint,related_files


def item(name):
    return {"id":"sampleid","attributes":{
        "kind":"file","name":name,"size":15,
        "extra":{"hashes":{"sha256":"a"*64,"md5":"b"*32}}
        },"links":{"download":"https://files.osf.io/v1/resources/a/b"}}


class OsfMetadataTests(unittest.TestCase):
    def test_rejects_redirect_or_non_osf_metadata_endpoint(self):
        for url in ("http://api.osf.io/v2/files/","https://evil.example/v2",
                    "https://api.osf.io.evil.example/v2"):
            with self.assertRaises(ValueError):endpoint(url)

    def test_related_folder_urls(self):
        r={"relationships":{"files":{"links":{"related":{
            "href":"https://api.osf.io/v2/files/abc/"}}}}}
        self.assertEqual(related_files(r),"https://api.osf.io/v2/files/abc/")

    def test_three_metadata_entries_without_eeg_reads(self):
        docs={"data":[item(f"sub-024_task-P3_{suffix}") for suffix
                      in ("eeg.set","eeg.fdt","events.tsv")],
              "links":{"next":None}}
        with patch("discover_holdout_osf_metadata.fetch_json",return_value=docs) as fake:
            r=discover("etdkz",[24])
        self.assertEqual(len(r["files"]),3)
        self.assertEqual(r["eeg_bytes_downloaded"],0)
        self.assertEqual(fake.call_count,1)

    def test_missing_file_fails_closed(self):
        with patch("discover_holdout_osf_metadata.fetch_json",
                   return_value={"data":[item("sub-024_task-P3_eeg.set")],
                                 "links":{"next":None}}):
            with self.assertRaisesRegex(ValueError,"missing"):
                discover("etdkz",[24])

    def test_duplicate_osf_name_fails_closed(self):
        d=item("sub-024_task-P3_eeg.set")
        with patch("discover_holdout_osf_metadata.fetch_json",
                   return_value={"data":[d,d],"links":{"next":None}}):
            with self.assertRaisesRegex(ValueError,"duplicate"):
                discover("etdkz",[24])


if __name__=="__main__":
    unittest.main()
