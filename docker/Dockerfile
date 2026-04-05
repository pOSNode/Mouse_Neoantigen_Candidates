FROM python:3.11-slim

LABEL maintainer="Owais Siddiqi <owais@owaissiddiqi.co.uk>"
LABEL description="RNA-Seq Neoantigen Discovery Pipeline"
LABEL version="1.0.0"

WORKDIR /app

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        tabix \
        bgzip \
        && rm -rf /var/lib/apt/lists/*

# Python dependencies
COPY docker/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy scripts
COPY scripts/ ./scripts/

# Non-root user for security
RUN useradd --create-home --shell /bin/bash pipeline
USER pipeline

CMD ["/bin/bash"]
