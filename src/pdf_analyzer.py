
from PyPDF2 import PdfReader

def analyze_pdf(pdf_path):
    """Analyze a PDF file and return its content and metadata."""
    reader = PdfReader(pdf_path)
    
    # Get basic metadata
    metadata = reader.metadata
    
    # Extract text from all pages
    text = ""
    for page in reader.pages:
        text += page.extract_text()
        
    return {
        "metadata": metadata,
        "text": text,
        "num_pages": len(reader.pages)
    }

if __name__ == "__main__":
    pdf_path = "src/docs/gdmpidc.pdf"
    analysis = analyze_pdf(pdf_path)
    print(f"Number of pages: {analysis['num_pages']}")
    print("\nMetadata:")
    for key, value in analysis['metadata'].items():
        print(f"{key}: {value}")
    print("\nFirst 500 characters of text:")
    print(analysis['text'][:500])
