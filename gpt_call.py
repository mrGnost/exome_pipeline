import argparse
import configparser
from pathlib import Path
from openai import OpenAI

PROMPT = (
    "You are a senior medical data analyst with high skills in bioinformatics and genomics. "
    "Look through the exomic file {vcf_filename}, which is annotated with clinvar, "
    "dbsnp and alphamissense, and decide which type of genetic disease the patient could have. "
    "Write only genes, possible diseases and a short supplementary info to the file gpt_report.txt."
)


def send_query(vcf_path: Path, client: OpenAI):
    prompt = PROMPT.format(vcf_filename=vcf_path.name)

    with open(vcf_path, "rb") as f:
        uploaded_file = client.files.create(file=f, purpose="assistants")

    response = client.chat.completions.create(
        model="gpt-5.4-2026-03-05",
        messages=[
            {
                "role": "user",
                "content": [
                    {"type": "text", "text": prompt},
                    {"type": "file", "file": {"file_id": uploaded_file.id}},
                ],
            }
        ]
    )

    Path("gpt_response.txt").write_text(response.choices[0].message.content)
    client.files.delete(uploaded_file.id)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_file", help="Path to the VCF file")
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read("config.ini")
    api_key = config["openai"]["api_key"]

    openai_client = OpenAI(api_key=api_key)
    send_query(Path(args.vcf_file), openai_client)
